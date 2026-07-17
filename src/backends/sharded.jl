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

function _build_rhs(
        eqs::MeanfieldEquations, ps, backend::ShardedBackend; jac::Union{Bool, Symbol} = false
    )
    mode = _sharded_jac_mode(jac)
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
    # bind through locals so `chunk`/`parallel` narrow to concrete `Int`/`Bool` (the union
    # field types otherwise flow into the `Int` chunk size and the `Bool` ShardedRHS slot)
    ch = backend.chunk
    chunk = ch isa Integer ? ch : max(1, cld(neq, 3 * nth))
    pl = backend.parallel
    parallel = pl isa Bool ? pl : (nth > 1 && neq >= SHARD_PARALLEL_MIN)
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
    rhs = ShardedRHS(fws, ranges, pars, pvec, parallel)
    mode === :none && return rhs
    return _attach_jacobian(rhs, mode, rhss, vars, codegen_pars, eqs.iv, neq)
end

# ---- Jacobian & time-gradient for stiff/implicit solvers (issue #294) -----------------
# The integrator disables its OWN sparse/colored finite differencing for complex states
# (its default engine is ForwardDiff, which has no complex support), so it would fall back
# to a dense O(n) column sweep. We instead hand it a ready sparse `jac` callable. Both modes
# compute the real-directional linearization (the derivative along the real perturbation
# axis, which is exactly what the integrator's own dense FD computes); this matches for
# conjugate-folded AND unfolded closures alike, so unlike the kernel's holomorphic M·v
# Jacobian there is no `HolomorphicJacobianError` here. A finite-difference `tgrad` is
# attached too, so Rosenbrock methods run without the caller passing
# `autodiff = AutoFiniteDiff()` to route around the same ForwardDiff-on-complex gap.

struct ShardedRHSWithJac{J, T}
    rhs::ShardedRHS
    jac::J                                       # (J, u, p, t) -> nothing, fills J.nzval
    tgrad::T                                     # (dT, u, p, t) -> nothing
    jac_prototype::SparseMatrixCSC{ComplexF64, Int}
end

# the ODEFunction carries the plain `ShardedRHS` as its `f` (so `prob.p`/`update_parameters!`
# keep working unchanged); the Jacobian and time-gradient ride alongside.
function _ode_function(w::ShardedRHSWithJac)
    return SciMLBase.ODEFunction{true}(
        w.rhs; jac = w.jac, tgrad = w.tgrad, jac_prototype = w.jac_prototype,
    )
end

function _sharded_jac_mode(jac)
    jac === false && return :none
    (jac === true || jac === :fd) && return :fd
    jac === :analytic && return :analytic
    return throw(
        ArgumentError(
            "ShardedBackend `jac` must be `false`, `true`/`:fd` (colored finite differences) " *
                "or `:analytic` (symbolic); got $(repr(jac)).",
        )
    )
end

function _attach_jacobian(rhs::ShardedRHS, mode, rhss, vars, codegen_pars, iv, neq)
    tgrad = _make_sharded_tgrad(rhs, neq)
    if mode === :analytic
        built = try
            _make_sharded_analytic_jac(rhss, vars, codegen_pars, iv, neq)
        catch err
            @info "ShardedBackend: analytic Jacobian construction failed; falling back to " *
                "colored finite differences." exception = err maxlog = 1
            nothing
        end
        built === nothing || return ShardedRHSWithJac(rhs, built[1], tgrad, built[2])
    end
    sp = _zeros_pattern(Symbolics.jacobian_sparsity(rhss, vars))
    colors = _greedy_column_colors(sp)
    return ShardedRHSWithJac(rhs, _make_sharded_fd_jac(rhs, sp, colors), tgrad, sp)
end

"""ComplexF64 sparse matrix carrying only `A`'s structure (all stored values zero)."""
_zeros_pattern(A::SparseMatrixCSC) =
    SparseMatrixCSC(size(A, 1), size(A, 2), copy(A.colptr), copy(A.rowval), zeros(ComplexF64, nnz(A)))

"""
Greedy distance-1 column coloring: two columns sharing a nonzero row get different colors,
so one finite-difference perturbation per color recovers every column of that color at once
(`ncolors + 1` RHS evaluations per Jacobian instead of `n + 1`).
"""
function _greedy_column_colors(A::SparseMatrixCSC)
    rows = rowvals(A)
    colors = zeros(Int, size(A, 2))
    rowcolors = [Set{Int}() for _ in axes(A, 1)]
    for j in axes(A, 2)
        used = Set{Int}()
        for k in nzrange(A, j)
            union!(used, rowcolors[rows[k]])
        end
        c = 1
        while c in used
            c += 1
        end
        colors[j] = c
        for k in nzrange(A, j)
            push!(rowcolors[rows[k]], c)
        end
    end
    return colors
end

# forward-difference step: sqrt(eps) scaled by the magnitude of the perturbed entry.
_fd_step(x) = sqrt(eps(Float64)) * max(abs(x), one(Float64))

"""
Colored forward-difference Jacobian filling `J`'s stored values (`J` shares `sp`'s CSC
pattern). Scratch is allocated per call, so the closure is reentrant on the same footing as
the RHS (safe to call concurrently, e.g. one Jacobian per trajectory).
"""
function _make_sharded_fd_jac(
        rhs::ShardedRHS, sp::SparseMatrixCSC{ComplexF64, Int}, colors::Vector{Int}
    )
    n = size(sp, 1)
    rows = rowvals(sp)
    groups = [findall(==(c), colors) for c in 1:maximum(colors; init = 0)]
    return function (J, u, p, t)
        f0 = Vector{ComplexF64}(undef, n)
        fp = Vector{ComplexF64}(undef, n)
        up = Vector{ComplexF64}(undef, n)
        rhs(f0, u, p, t)
        nzv = nonzeros(J)
        @inbounds for g in groups
            copyto!(up, u)
            for j in g
                up[j] += _fd_step(u[j])
            end
            rhs(fp, up, p, t)
            for j in g
                h = _fd_step(u[j])
                for k in nzrange(sp, j)
                    nzv[k] = (fp[rows[k]] - f0[rows[k]]) / h
                end
            end
        end
        return nothing
    end
end

"""Forward-difference time-gradient ∂f/∂t (identically zero for autonomous drifts)."""
function _make_sharded_tgrad(rhs::ShardedRHS, n::Int)
    return function (dT, u, p, t)
        f0 = Vector{ComplexF64}(undef, n)
        fp = Vector{ComplexF64}(undef, n)
        rhs(f0, u, p, t)
        h = _fd_step(t)
        rhs(fp, u, p, t + h)
        @inbounds for i in 1:n
            dT[i] = (fp[i] - f0[i]) / h
        end
        return nothing
    end
end

"""
Real-directional analytic Jacobian of `rhss` wrt `vars`. The sharded drifts express a
moment's conjugate partner as `conj(vₖ)`, which is not holomorphic, so `Symbolics.jacobian`
would leave `d conj(vₖ)/d vₖ` unevaluated. Replacing each `conj(vₖ)` by a fresh `wₖ`,
differentiating wrt both `vₖ` and `wₖ` and summing the columns gives the derivative along
the real perturbation axis (`d conj(v)/d(Re v) = 1`), i.e. the same linearization the
integrator's finite differencing computes; then `wₖ` is restored to `conj(vₖ)`.
"""
function _real_directional_jacobian(rhss, vars)
    w = [Symbolics.unwrap(Symbolics.variable(Symbol(:__qc_conj, i))) for i in eachindex(vars)]
    fwd = Dict{Any, Any}(Symbolics.unwrap(conj(vars[i])) => w[i] for i in eachindex(vars))
    rsub = Any[Symbolics.substitute(r, fwd) for r in rhss]
    J = Symbolics.jacobian(rsub, vars) .+ Symbolics.jacobian(rsub, w)
    back = Dict{Any, Any}(w[i] => Symbolics.unwrap(conj(vars[i])) for i in eachindex(vars))
    return map(e -> Symbolics.substitute(e, back), J)
end

"""
Analytic sparse Jacobian: build the real-directional Jacobian symbolically, drop structural
zeros, and codegen an in-place filler with the same `(u, p, t)` signature as the RHS chunks.
Returns `(jac!, jac_prototype)`.
"""
function _make_sharded_analytic_jac(rhss, vars, codegen_pars, iv, n)
    Jsym = _real_directional_jacobian(rhss, vars)
    is = Int[]
    js = Int[]
    vs = Symbolics.Num[]
    for j in 1:n, i in 1:n
        e = Jsym[i, j]
        SymbolicUtils._iszero(Symbolics.unwrap(e)) && continue
        push!(is, i)
        push!(js, j)
        push!(vs, Symbolics.Num(e))
    end
    Jsp = sparse(is, js, vs, n, n)
    expr = Symbolics.build_function(Jsp, vars, codegen_pars, iv; expression = Val{true}, cse = true)[2]
    jf = @RuntimeGeneratedFunction(expr)
    proto = _zeros_pattern(Jsp)
    jac = (J, u, p, t) -> (jf(J, u, p, t); nothing)
    return jac, proto
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
