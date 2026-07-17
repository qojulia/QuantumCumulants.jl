# Sharded native codegen backend (issue #294): the completed drifts are compiled to
# chunked `build_function` kernels (RuntimeGeneratedFunctions), precompiled in parallel,
# and driven through a FunctionWrapper table so the driver loop is one concrete type.
# Handles anything `Symbolics.build_function` can compile (t-dependent and non-polynomial
# drifts included), at a higher cold-construction cost than the kernel backend.

const VT = SubArray{ComplexF64, 1, Vector{ComplexF64}, Tuple{UnitRange{Int64}}, true}
const FWT = FunctionWrapper{Nothing, Tuple{VT, Vector{ComplexF64}, Vector{ComplexF64}, Float64}}

struct ShardedRHS
    fws::Vector{FWT}
    ranges::Vector{UnitRange{Int}}
    params::Vector{Any}
    pvec::Vector{ComplexF64}     # the live parameter vector handed out as prob.p
    parallel::Bool               # run the chunks concurrently (disjoint du ranges)
end
function (r::ShardedRHS)(du::Vector{ComplexF64}, u::Vector{ComplexF64}, p, t)
    if r.parallel
        # chunks write disjoint `du` ranges and only read shared `u`/`p`, so the loop is
        # data-race free without locking; :dynamic load-balances the uneven chunk costs.
        Threads.@threads :dynamic for i in eachindex(r.fws)
            @inbounds r.fws[i](view(du, r.ranges[i]), u, p, t)
        end
    else
        @inbounds for i in eachindex(r.fws)
            r.fws[i](view(du, r.ranges[i]), u, p, t)
        end
    end
    return nothing
end
function (r::ShardedRHS)(du, u, p, t)
    throw(
        ArgumentError(
            "the sharded RHS is compiled for ComplexF64 states; got eltype $(eltype(u)). " *
                "Implicit solvers with ForwardDiff autodiff are unsupported here; use an explicit " *
                "solver or `KernelBackend()` with `jac = true`.",
        )
    )
end

"""
Drifts resolved directly off the graph (no `System`/`mtkcompile` detour), through the
SAME treatments-aware machinery as `MTK.System`.
"""
function _sharded_expressions(eqs)
    reg = _state_registry(eqs)
    unresolved = Set{Any}()
    resolve = _leaf_resolver(reg, unresolved)
    rhss = Any[
        SymbolicUtils.unwrap(mapleaves(resolve, SymbolicUtils.unwrap(eq.rhs)))
            for eq in eqs.equations
    ]
    isempty(unresolved) || throw(UnresolvedMomentError(first(unresolved)))
    return rhss, reg.vars
end

"""
Scalar parameter occurrences of the drifts: `discover_params` over the drift expressions
(t allowed here, so no iv check), minus the state variables themselves.
"""
function _sharded_params(rhss, vars, iv_uw)
    statelike = Set{Any}(vars)
    pars = Any[]
    for p in discover_params(rhss)
        (p in statelike || isequal(p, iv_uw)) && continue
        push!(pars, p)
    end
    return pars
end

function _build_rhs(eqs::MeanfieldEquations, ps, backend::ShardedBackend; jac::Bool = false)
    jac && throw(
        ArgumentError(
            "`jac = true` is not supported by ShardedBackend (colored-FD Jacobian not yet " *
                "prototyped); use KernelBackend.",
        )
    )
    iv_uw = SymbolicUtils.unwrap(eqs.iv)
    rhss, vars = _sharded_expressions(eqs)
    pars = _sharded_params(rhss, vars, iv_uw)
    vals = kernel_pdict(pars, parameter_map(eqs, ps))
    # array accesses (evaluated systems) are aliased to plain scalars for codegen
    nonsym = Any[p for p in pars if !SymbolicUtils.issym(p)]
    if !isempty(nonsym)
        aliases = Dict{Any, Any}(
            p => Symbolics.unwrap(Symbolics.variable(Symbol(:__qc_p, i)))
                for (i, p) in enumerate(nonsym)
        )
        rhss = Any[Symbolics.substitute(r, aliases) for r in rhss]
        codegen_pars = Any[get(aliases, p, p) for p in pars]
    else
        codegen_pars = pars
    end
    pvec = ComplexF64[ComplexF64(vals[p]) for p in pars]
    neq = length(rhss)
    nth = clamp(backend.threads, 1, Threads.nthreads())
    chunk = backend.chunk === :auto ? max(1, cld(neq, 3 * nth)) : backend.chunk
    parallel = backend.parallel === :auto ? (nth > 1 && neq >= SHARD_PARALLEL_MIN) :
        backend.parallel
    ranges = UnitRange{Int}[r[1]:r[end] for r in Iterators.partition(1:neq, chunk)]
    # `build_function` Expr generation is pure and heterogeneous across chunks; run it in
    # parallel (was serial) on `nth` tasks. `tmap` returns Exprs in `ranges` order.
    # `@RuntimeGeneratedFunction` construction stays serial (global RGF registration).
    # The `let` localizes `rhss`/`codegen_pars` (reassigned above, hence boxed) so the task
    # closure does not capture boxed variables, which OhMyThreads rejects.
    fexprs = let rhss = rhss, vars = vars, codegen_pars = codegen_pars, iv = eqs.iv
        tmap(
            r -> Symbolics.build_function(
                rhss[r], vars, codegen_pars, iv; expression = Val{true}, cse = true,
            )[2],
            Expr, ranges; scheduler = DynamicScheduler(; ntasks = nth),
        )
    end
    fs = [@RuntimeGeneratedFunction(fe) for fe in fexprs]
    sig = (VT, Vector{ComplexF64}, Vector{ComplexF64}, Float64)
    # warm up each chunk's native code; per-chunk compile cost varies widely, so a greedy
    # schedule balances better than the round-robin static split it replaces.
    tforeach(fs; scheduler = GreedyScheduler(; ntasks = nth)) do f
        precompile(f, sig)
    end
    fws = FWT[
        FWT(
                let g = f
                    (du, u, p, t) -> (g(du, u, p, t); nothing)
            end
            ) for f in fs
    ]
    return ShardedRHS(fws, ranges, pars, pvec, parallel)
end

_prob_p(r::ShardedRHS) = r.pvec

update_parameters!(r::ShardedRHS, pdict) = update_parameters!(r, r.pvec, pdict)
function update_parameters!(r::ShardedRHS, p::AbstractVector, pdict)
    fresh = kernel_pdict(r.params, pdict; strict = false)
    for (i, par) in enumerate(r.params)
        haskey(fresh, par) && (p[i] = fresh[par])
    end
    return r
end
