# Public ODE surface (issue #294): `ODEFunction`/`ODEProblem` constructors on
# `MeanfieldEquations` (QC owns the type, no piracy). The `System(eqs)`/MTK route
# stays untouched for interop; these constructors compile the RHS directly.

"""
    RHSBackend

Abstract supertype of the RHS compilation strategies accepted by
`ODEFunction(eqs, ps; backend)` and `ODEProblem(eqs, u0, tspan, ps; backend)`.
"""
abstract type RHSBackend end

"""
    KernelBackend(; cache = nothing, parallel = :auto)

Compile the RHS to the moment-polynomial kernel `du = M * v`: the completed drifts are
lowered once to sparse data (no per-model native code), so construction is fast and the
RHS allocation-free. Requires drifts polynomial in the moments (anything
`meanfield`/`complete!` produces). Pass `cache = path` to persist the lowered tables in
a JLD2 file keyed by a digest of the equations (requires `using JLD2`).

The sparse `M * v` is stored in CSR layout and evaluated as an independent per-row gather,
which alone is several times faster than a CSC `mul!`. `parallel = :auto` additionally
threads both RHS passes (the monomial update and the gather) through Polyester's persistent
task pool on systems of `≥ $(KERNEL_PARALLEL_MIN)` equations running on `> 1` threads, for a
2-3x further speedup; pass `true`/`false` to force it. Serial and threaded produce
bit-identical results (every product and row-sum runs in the same order).
"""
struct KernelBackend <: RHSBackend
    cache::Union{Nothing, String}
    parallel::Union{Bool, Symbol}
end
KernelBackend(; cache = nothing, parallel = :auto) = KernelBackend(cache, parallel)

# Equation count above which the parallel driver's ~4 µs threading overhead pays off
# (measured against the evaluated superradiant laser on 12 threads: break-even near
# ~400 eqs, a clear win by ~1000). Only consulted for `parallel = :auto`.
const SHARD_PARALLEL_MIN = 512

"""
    ShardedBackend(; chunk = :auto, threads = Threads.nthreads(), parallel = :auto)

Compile the RHS to native code in chunks of `chunk` equations (RuntimeGeneratedFunctions
behind a FunctionWrapper table), precompiled on `threads` threads. Handles anything
`Symbolics.build_function` can compile, including t-dependent and non-polynomial drifts;
cold construction is slower than `KernelBackend` and grows with system size.

`chunk = :auto` (the default) sizes the chunks to `≈ 3·threads` chunks total, which both
balances the parallel precompile and keeps the runtime driver's per-chunk overhead low; an
integer forces a fixed chunk size. `parallel = :auto` runs the chunks of the compiled RHS
concurrently (`Threads.@threads :dynamic`, chunks write disjoint `du` ranges) once the
system is large enough that the ~4 µs threading overhead pays off (`≥ $(SHARD_PARALLEL_MIN)`
equations on `> 1` threads); pass `true`/`false` to force it either way.
"""
struct ShardedBackend <: RHSBackend
    chunk::Union{Int, Symbol}
    threads::Int
    parallel::Union{Bool, Symbol}
end
ShardedBackend(; chunk = :auto, threads = Threads.nthreads(), parallel = :auto) =
    ShardedBackend(chunk, threads, parallel)

"""
    AutoBackend()

The default: try `KernelBackend()`, falling back to `ShardedBackend()` exactly when the
drift is not polynomial in the moments (`NonPolynomialDriftError`). Every other lowering
error stays a hard error.
"""
struct AutoBackend <: RHSBackend end

struct KernelRHS{F, J}
    kernel::MomentKernel
    kp::KernelParameters{F}
    jacobian::J                  # `nothing`, or a `JacKernel` sharing the kernel's `v`
end
KernelRHS(kernel::MomentKernel, kp::KernelParameters) = KernelRHS(kernel, kp, nothing)
(r::KernelRHS)(du, u::AbstractVector{ComplexF64}, p::KernelParameters, t) =
    r.kernel(du, u, p, t)
function (r::KernelRHS)(du, u, p::KernelParameters, t)
    throw(
        ArgumentError(
            "the moment kernel is compiled for ComplexF64 states; got eltype $(eltype(u)). " *
                "Implicit solvers with ForwardDiff autodiff are unsupported on this RHS; use an " *
                "explicit solver, or `jac = true` for the analytic sparse Jacobian.",
        )
    )
end
function (r::KernelRHS)(du, u, p, t)
    throw(
        ArgumentError(
            "a kernel ODEProblem carries a `KernelParameters` object as `prob.p`; got " *
                "$(typeof(p)). Use `update_parameters!(prob, Dict(...))` instead of " *
                "`remake(prob, p = ...)`.",
        )
    )
end

"""
Normalize the `jac` keyword for the kernel backend to a `Bool`. The kernel's Jacobian is
always analytic (the M·v structure), so `true`/`:analytic` request it and `:fd` is rejected
(finite differences are a `ShardedBackend` option).
"""
function _kernel_jac_flag(jac)
    jac === false && return false
    (jac === true || jac === :analytic) && return true
    return throw(
        ArgumentError(
            "KernelBackend's Jacobian is analytic; `jac` must be `false`, `true` or " *
                "`:analytic`, got $(repr(jac)). Use `ShardedBackend` for `:fd`.",
        )
    )
end

function _build_rhs(
        eqs::MeanfieldEquations, ps, backend::KernelBackend; jac::Union{Bool, Symbol} = false
    )
    dojac = _kernel_jac_flag(jac)
    # cache flow (jac = true relowers fresh: the stored tables lack the Jacobian's
    # complement monomials, so a hit would be structurally incomplete)
    text = digest = nothing
    if backend.cache !== nothing && !dojac
        text = canonical_text(eqs)
        digest = bytes2hex(sha256(text))
        entry = _cache_load(backend.cache, text, digest)
        entry === nothing || return _loaded_rhs(entry, eqs, ps, backend.parallel)
    end
    ir = lower(eqs)
    backend.cache === nothing || dojac || _cache_store!(backend.cache, text, digest, ir)
    values = kernel_pdict(ir.params, parameter_map(eqs, ps))
    cvals = coefficient_values(ir, values)
    par = _resolve_kernel_parallel(backend.parallel, ir.nstates)
    if !dojac
        kernel = MomentKernel(ir, cvals; parallel = par)
        return KernelRHS(kernel, KernelParameters(ir, kernel.Mt, values))
    end
    # extended IR so `v` covers the delete-one complement monomials; RHS and Jacobian
    # share one monomial id space (and the same workspace `v`)
    ir_ext, jir = jacobian_ir(ir)
    kernel = MomentKernel(ir_ext, cvals; parallel = par)
    jk = JacKernel(jir, kernel.parent, kernel.leaf, cvals, kernel.v)
    return KernelRHS(kernel, KernelParameters(ir_ext, kernel.Mt, values), jk)
end

"""
Kernel RHS from a cache entry: plain tables plus the stored coefficient evaluator RGF.
Parameters match by printed name (`kernel_pdict` on the stored name strings), so a
loaded kernel is fully sweepable through the same `KernelParameters` machinery.
"""
function _loaded_rhs(entry, eqs, ps, parallel_flag)
    values = kernel_pdict(Any[entry.params...], parameter_map(eqs, ps))
    evalcoeffs = vals -> ComplexF64.(entry.rgf(ComplexF64[vals[n] for n in entry.params]))
    cvals = evalcoeffs(values)
    Mt = sparse(entry.coo_j, entry.coo_i, cvals[entry.coo_c], length(entry.parent), entry.nstates, +)
    v = _make_vbufs(length(entry.parent))
    par = _resolve_kernel_parallel(parallel_flag, entry.nstates)
    kernel = MomentKernel(Mt, entry.parent, entry.leaf, v, par)
    kp = KernelParameters(
        Any[entry.params...], values, evalcoeffs, entry.coo_c,
        nz_map(Mt, entry.coo_j, entry.coo_i),
    )
    return KernelRHS(kernel, kp)
end

function _build_rhs(eqs::MeanfieldEquations, ps, ::AutoBackend; jac::Union{Bool, Symbol} = false)
    return try
        # the kernel's only Jacobian is analytic, so any request maps to `true` there
        _build_rhs(eqs, ps, KernelBackend(); jac = (jac === false ? false : true))
    catch e
        e isa NonPolynomialDriftError || rethrow()
        @info "drift is not polynomial in the moments; falling back to ShardedBackend " *
            "(chunked native codegen). Cold construction will be slower; `jac` uses finite " *
            "differences unless `:analytic` is requested." maxlog = 1
        _build_rhs(eqs, ps, ShardedBackend(); jac)
    end
end

"""
    SciMLBase.ODEFunction(eqs::MeanfieldEquations, ps; backend = AutoBackend(), jac = false)

In-place `ODEFunction` over the completed moment equations, compiled by `backend`.
Parameters are required at construction; `ps` is a `Dict` or pairs of `parameter => value`.
With `jac = true` the function carries a sparse Jacobian and its `jac_prototype`: the kernel
backend attaches its analytic M·v Jacobian, while `ShardedBackend` attaches a colored
finite-difference Jacobian (`jac = :analytic` there requests a symbolic one instead) plus a
finite-difference time-gradient, so implicit solvers run on the complex state directly.
"""
function SciMLBase.ODEFunction(
        eqs::MeanfieldEquations, ps; backend::RHSBackend = AutoBackend(),
        jac::Union{Bool, Symbol} = false,
    )
    rhs = _build_rhs(eqs, ps, backend; jac)
    return _ode_function(rhs)
end

_ode_function(rhs) = SciMLBase.ODEFunction{true}(rhs)
function _ode_function(rhs::KernelRHS)
    rhs.jacobian === nothing && return SciMLBase.ODEFunction{true}(rhs)
    return SciMLBase.ODEFunction{true}(
        rhs; jac = rhs.jacobian, jac_prototype = copy(rhs.jacobian.jac.Jproto),
    )
end

_prob_p(r::KernelRHS) = r.kp

"""
    SciMLBase.ODEProblem(eqs::MeanfieldEquations, u0, tspan, ps; backend = AutoBackend(), kwargs...)

Build the `ODEProblem` directly from the completed moment equations. `u0` is a numeric
vector aligned with `eqs.states`, an `AbstractDict` keyed by state averages, or a numeric
quantum state (ket / density operator). `ps` maps every symbolic parameter to its value.
"""
function SciMLBase.ODEProblem(
        eqs::MeanfieldEquations, u0, tspan, ps;
        backend::RHSBackend = AutoBackend(), jac::Union{Bool, Symbol} = false, kwargs...,
    )
    f = SciMLBase.ODEFunction(eqs, ps; backend, jac)
    return SciMLBase.ODEProblem(f, _u0_vector(eqs, u0), tspan, _prob_p(f.f); kwargs...)
end

function _u0_vector(eqs, u0::AbstractVector{<:Number})
    length(u0) == length(eqs.states) || throw(
        DimensionMismatch(
            "u0 length $(length(u0)) does not match number of states $(length(eqs.states))",
        )
    )
    return ComplexF64.(u0)
end
function _u0_vector(eqs, u0::AbstractDict)
    reg = _state_registry(eqs)
    return ComplexF64[
        haskey(u0, reg.vars[k]) ? ComplexF64(u0[reg.vars[k]]) :
            ComplexF64(get(u0, eqs.states[k], 0)) for k in eachindex(eqs.states)
    ]
end
_u0_vector(eqs, state) = initial_values(eqs, state)   # kets / density operators

# ---- parameter sweeps and ensembles ---------------------------------------------------

"""
    update_parameters!(prob, pdict)
    update_parameters!(f::ODEFunction, pdict)

Rewrite the compiled RHS in place for new parameter values (a `Dict` or pairs of
`parameter => value`; unmentioned parameters keep their current values). This is the
sweep path: the lowered tables are reused, only the numeric coefficients are refreshed.
"""
function update_parameters!(r::KernelRHS, pdict)
    kp = r.kp
    # resolve the user's dict on its own (arrays expand to their per-slot accesses), so a
    # fresh value always overrides the stored one; undetermined parameters keep theirs
    fresh = kernel_pdict(kp.params, pdict; strict = false)
    merge!(kp.values, fresh)
    cvals = kp.evalcoeffs(kp.values)
    write_nzval!(r.kernel, kp, cvals)
    r.jacobian === nothing || copyto!(r.jacobian.c, cvals)
    return r
end

function Base.copy(r::KernelRHS)
    k = r.kernel
    # structure tables (parent/leaf/fac/fac_ptr) are read-only and shared; own Mt and v
    # (deep-copy the per-thread scratch set so the copy's buffers are independent)
    k2 = MomentKernel(copy(k.Mt), k.parent, k.leaf, k.fac, k.fac_ptr, map(copy, k.v), k.parallel)
    kp = r.kp
    kp2 = KernelParameters(
        kp.params, Dict{Any, Any}(kp.values), kp.evalcoeffs, kp.coo_c, kp.nzmap
    )
    r.jacobian === nothing && return KernelRHS(k2, kp2)
    # own coefficient values, shared structure tables, workspace shared with the copy's kernel
    jk = r.jacobian
    return KernelRHS(k2, kp2, JacKernel(jk.jac, k2.parent, k2.leaf, copy(jk.c), k2.v))
end

function update_parameters!(prob::SciMLBase.ODEProblem, pdict)
    r = prob.f.f
    if r isa ShardedRHS
        # write the vector the problem actually carries (a remade prob.p copy included)
        update_parameters!(r, prob.p, pdict)
    else
        update_parameters!(prob.f, pdict)
    end
    return prob
end
update_parameters!(f::SciMLBase.ODEFunction, pdict) = (update_parameters!(f.f, pdict); f)

function Base.copy(f::SciMLBase.ODEFunction{iip, spec, <:KernelRHS}) where {iip, spec}
    r2 = copy(f.f)
    r2.jacobian === nothing && return SciMLBase.ODEFunction{iip}(r2)
    return SciMLBase.ODEFunction{iip}(
        r2; jac = r2.jacobian, jac_prototype = copy(f.jac_prototype),
    )
end

# ---- guard methods: typed errors on QC types, no fall-through to SciMLBase -------------

SciMLBase.ODEProblem(eqs::AbstractMeanfieldEquations, u0, tspan) = throw(
    ArgumentError(
        "parameters are required at construction: `ODEProblem(eqs, u0, tspan, ps; ...)` with " *
            "`ps` a Dict or pairs of parameter => value (the kernel backend bakes them into " *
            "its tables).",
    )
)
SciMLBase.ODEFunction(eqs::AbstractMeanfieldEquations) = throw(
    ArgumentError("parameters are required at construction: `ODEFunction(eqs, ps; ...)`.")
)
SciMLBase.ODEProblem(eqs::NoiseMeanfieldEquations, u0, tspan, ps; kwargs...) = throw(
    ArgumentError(
        "noise systems are SDEs; the fast RHS backends are deterministic-only in v1. Build " *
            "the MTK path instead: `System(eqs; name = ...)`.",
    )
)
