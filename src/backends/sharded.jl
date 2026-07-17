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
end
function (r::ShardedRHS)(du::Vector{ComplexF64}, u::Vector{ComplexF64}, p, t)
    @inbounds for i in eachindex(r.fws)
        r.fws[i](view(du, r.ranges[i]), u, p, t)
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
    ranges = UnitRange{Int}[r[1]:r[end] for r in Iterators.partition(1:neq, backend.chunk)]
    fexprs = [
        Symbolics.build_function(
                rhss[r], vars, codegen_pars, eqs.iv; expression = Val{true}, cse = true,
            )[2] for r in ranges
    ]
    fs = [@RuntimeGeneratedFunction(fe) for fe in fexprs]
    sig = (VT, Vector{ComplexF64}, Vector{ComplexF64}, Float64)
    nth = clamp(backend.threads, 1, Threads.nthreads())
    tasks = [
        Threads.@spawn(
                for f in fs[i:nth:end]
                    precompile(f, sig)
            end
            ) for i in 1:nth
    ]
    foreach(wait, tasks)
    fws = FWT[
        FWT(
                let g = f
                    (du, u, p, t) -> (g(du, u, p, t); nothing)
            end
            ) for f in fs
    ]
    return ShardedRHS(fws, ranges, pars, pvec)
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
