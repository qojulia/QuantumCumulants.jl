# Direct SciML surface for the compact structured numerical evaluator.

"""
    RHSBackend

Abstract supertype of direct RHS execution strategies. The first production strategy is
`KernelBackend`; unsupported systems should use `System(eqs)` and the ModelingToolkit path.
"""
abstract type RHSBackend end

"""
    KernelBackend(; parallel = :auto)

Lower completed `MeanfieldEquations` to the compact moment-polynomial representation
``du = M * v``. The lowering is structural and is independent of the numeric parameter
values. `parallel = :auto` enables the Polyester evaluation path for sufficiently large
systems on a multi-threaded Julia process; pass `true` or `false` to select it explicitly.

The direct evaluator currently uses `ComplexF64` state and coefficient tables.
"""
struct KernelBackend <: RHSBackend
    parallel::Union{Bool,Symbol}
end

function KernelBackend(; parallel = :auto)
    parallel === :auto ||
        parallel isa Bool ||
        throw(ArgumentError("KernelBackend(parallel) must be `:auto`, `true`, or `false`."))
    return KernelBackend(parallel)
end

struct KernelRHS{F,J}
    kernel::MomentKernel
    kp::KernelParameters{F}
    jacobian::J                  # `nothing`, or the analytic `JacKernel`
end
KernelRHS(kernel::MomentKernel, kp::KernelParameters) = KernelRHS(kernel, kp, nothing)

(r::KernelRHS)(du, u::AbstractVector{ComplexF64}, p::KernelParameters, t) =
    r.kernel(du, u, p, t)
function (r::KernelRHS)(du, u, p::KernelParameters, t)
    throw(
        ArgumentError(
            "the direct moment evaluator is compiled for ComplexF64 states; got " *
            "eltype $(eltype(u)). Use an explicit solver or provide `jac = true` " *
            "for the analytic sparse Jacobian.",
        ),
    )
end
function (r::KernelRHS)(du, u, p, t)
    throw(
        ArgumentError(
            "a direct ODEProblem carries a `KernelParameters` object as `prob.p`; got " *
            "$(typeof(p)). Use `update_parameters!(prob, Dict(...))` for sweeps.",
        ),
    )
end

"""Normalize the strict analytic Jacobian mode of the compact evaluator."""
function _kernel_jac_flag(jac)
    jac === false && return false
    (jac === true || jac === :analytic) && return true
    return throw(
        ArgumentError(
            "KernelBackend's Jacobian is analytic; `jac` must be `false`, `true`, or " *
            "`:analytic`, got $(repr(jac)). The direct backend does not provide FD " *
            "Jacobians; use `System(eqs)` for that path.",
        ),
    )
end

function _build_rhs(
    eqs::MeanfieldEquations,
    ps,
    backend::KernelBackend;
    jac::Union{Bool,Symbol} = false,
)
    dojac = _kernel_jac_flag(jac)
    ir = _lower_moment_ir(eqs)
    values = kernel_pdict(ir.params, parameter_map(eqs, ps))
    cvals = coefficient_values(ir, values)
    parallel = _resolve_kernel_parallel(backend.parallel, ir.nstates)
    if !dojac
        kernel = MomentKernel(ir, cvals; parallel)
        return KernelRHS(kernel, KernelParameters(ir, kernel.Mt, values))
    end

    # Extend the shared monomial table with delete-one complements needed by the analytic
    # Jacobian. `jacobian_ir` rejects any conjugate factor before this evaluator is built.
    ir_ext, jir = jacobian_ir(ir)
    kernel = MomentKernel(ir_ext, cvals; parallel)
    jacobian = JacKernel(ir_ext, jir, cvals)
    return KernelRHS(kernel, KernelParameters(ir_ext, kernel.Mt, values), jacobian)
end

"""
    SciMLBase.ODEFunction(eqs::MeanfieldEquations, ps; backend = KernelBackend(), jac = false)

Build an in-place ODE function directly from completed moment equations. `ps` is required at
construction because it initializes the numeric coefficient table. Use `jac = true` or
`:analytic` only for holomorphic closures.
"""
function SciMLBase.ODEFunction(
    eqs::MeanfieldEquations,
    ps;
    backend::RHSBackend = KernelBackend(),
    jac::Union{Bool,Symbol} = false,
)
    rhs = _build_rhs(eqs, ps, backend; jac)
    return _ode_function(rhs)
end

_ode_function(rhs) = SciMLBase.ODEFunction{true}(rhs)
function _ode_function(rhs::KernelRHS)
    rhs.jacobian === nothing && return SciMLBase.ODEFunction{true}(rhs)
    return SciMLBase.ODEFunction{true}(
        rhs;
        jac = rhs.jacobian,
        jac_prototype = copy(rhs.jacobian.jac.Jproto),
    )
end

_prob_p(r::KernelRHS) = r.kp

"""
    SciMLBase.ODEProblem(eqs::MeanfieldEquations, u0, tspan, ps; kwargs...)

Build an `ODEProblem` directly from completed equations. `u0` may be a `ComplexF64`-compatible
vector aligned with `eqs.states`, a state-keyed dictionary, or a numeric quantum state.
"""
function SciMLBase.ODEProblem(
    eqs::MeanfieldEquations,
    u0,
    tspan,
    ps;
    backend::RHSBackend = KernelBackend(),
    jac::Union{Bool,Symbol} = false,
    kwargs...,
)
    f = SciMLBase.ODEFunction(eqs, ps; backend, jac)
    return SciMLBase.ODEProblem(f, _u0_vector(eqs, u0), tspan, _prob_p(f.f); kwargs...)
end

function _u0_vector(eqs, u0::AbstractVector{<:Number})
    length(u0) == length(eqs.states) || throw(
        DimensionMismatch(
            "u0 length $(length(u0)) does not match number of states $(length(eqs.states))",
        ),
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
_u0_vector(eqs, state) = initial_values(eqs, state)

# ---- parameter sweeps and copies -------------------------------------------------------

"""
    update_parameters!(prob, pdict)
    update_parameters!(f::ODEFunction, pdict)

Refresh only the numeric coefficient data of a direct evaluator. The structural lowering and
monomial table are reused, so repeated solves and parameter sweeps do not rederive equations.
"""
function update_parameters!(r::KernelRHS, pdict)
    kp = r.kp
    fresh = kernel_pdict(kp.params, pdict; strict = false)
    merge!(kp.values, fresh)
    cvals = kp.evalcoeffs(kp.values)
    write_nzval!(r.kernel, kp, cvals)
    r.jacobian === nothing || copyto!(r.jacobian.c, cvals)
    return r
end

function Base.copy(r::KernelRHS)
    k = r.kernel
    # Structure tables are immutable after construction. Own the sparse values and every
    # scratch buffer so copied problems can be updated and evaluated independently.
    k2 = MomentKernel(copy(k.Mt), k.parent, k.leaf, k.fac, k.fac_ptr, copy(k.v), k.parallel)
    kp = r.kp
    kp2 = KernelParameters(
        kp.params,
        Dict{Any,Any}(kp.values),
        kp.evalcoeffs,
        kp.coo_c,
        kp.nzmap,
    )
    r.jacobian === nothing && return KernelRHS(k2, kp2)
    jk = r.jacobian
    jac2 = JacKernel(jk.jac, k2.parent, k2.leaf, copy(jk.c), copy(jk.v))
    return KernelRHS(k2, kp2, jac2)
end

function update_parameters!(prob::SciMLBase.ODEProblem, pdict)
    r = prob.f.f
    r isa KernelRHS || throw(
        ArgumentError(
            "update_parameters! is only available for direct KernelBackend problems.",
        ),
    )
    update_parameters!(r, pdict)
    return prob
end
update_parameters!(f::SciMLBase.ODEFunction, pdict) = (update_parameters!(f.f, pdict); f)

function Base.copy(f::SciMLBase.ODEFunction{iip,spec,<:KernelRHS}) where {iip,spec}
    r2 = copy(f.f)
    r2.jacobian === nothing && return SciMLBase.ODEFunction{iip}(r2)
    return SciMLBase.ODEFunction{iip}(
        r2;
        jac = r2.jacobian,
        jac_prototype = copy(f.jac_prototype),
    )
end

# ---- guard methods ---------------------------------------------------------------------

SciMLBase.ODEProblem(eqs::AbstractMeanfieldEquations, u0, tspan) = throw(
    ArgumentError(
        "parameters are required at construction: `ODEProblem(eqs, u0, tspan, ps; ...)` " *
        "with `ps` a Dict or collection of parameter => value pairs.",
    ),
)
SciMLBase.ODEFunction(eqs::AbstractMeanfieldEquations) = throw(
    ArgumentError("parameters are required at construction: `ODEFunction(eqs, ps; ...)`."),
)
SciMLBase.ODEFunction(eqs::NoiseMeanfieldEquations, ps; kwargs...) = throw(
    ArgumentError(
        "noise systems are SDEs; direct ODE execution supports deterministic " *
        "`MeanfieldEquations` only. Build `System(eqs; name = ...)` for the MTK path.",
    ),
)
SciMLBase.ODEProblem(eqs::NoiseMeanfieldEquations, u0, tspan, ps; kwargs...) = throw(
    ArgumentError(
        "noise systems are SDEs; direct ODE execution supports deterministic " *
        "`MeanfieldEquations` only. Build `System(eqs; name = ...)` for the MTK path.",
    ),
)
