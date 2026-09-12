# SciML integration and cold parameter binding for the compact MomentKernel.
#
# The numerical boundary is deliberate:
#   * KernelParameterPlan retains symbolic construction metadata.
#   * KernelParameters is the concrete numeric `prob.p` payload.
#   * KernelRHS touches only `p.coeffs` on the hot path.

"""
    KernelBackend()

Explicit opt-in backend for constructing an `ODEProblem` directly from completed deterministic
`MeanfieldEquations` through `MomentIR` and the compact serial `MomentKernel`.
"""
struct KernelBackend end

struct KernelParameterPlan
    params::Vector{Any}
    coeffs::Vector{Any}
end

KernelParameterPlan(ir::MomentIR) = KernelParameterPlan(ir.params, ir.coeffs)

"""Concrete numerical parameter state carried as `prob.p` by direct kernel problems."""
struct KernelParameters{T}
    values::Vector{T}
    coeffs::Vector{T}
end

Base.copy(p::KernelParameters{T}) where {T} =
    KernelParameters{T}(copy(p.values), copy(p.coeffs))

struct KernelRHS{K, P, J}
    kernel::K
    plan::P
    jacobian::J
end

KernelRHS(kernel, plan) = KernelRHS(kernel, plan, nothing)

function (rhs::KernelRHS)(du, u, p::KernelParameters, t)
    rhs.kernel(du, u, p.coeffs)
    return nothing
end

# ---- parameter binding ---------------------------------------------------------------

function _kernel_param_name(p)
    if p isa SymbolicUtils.BasicSymbolic &&
            SymbolicUtils.iscall(p) &&
            SymbolicUtils.operation(p) === getindex
        args = SymbolicUtils.arguments(p)
        isempty(args) && return nothing
        base = args[1]
        return base isa SymbolicUtils.BasicSymbolic ? Base.nameof(base) : nothing
    end
    return _param_name(p)
end

function _kernel_param_slots(p)
    if p isa SymbolicUtils.BasicSymbolic &&
            SymbolicUtils.iscall(p) &&
            SymbolicUtils.operation(p) === getindex
        args = SymbolicUtils.arguments(p)
        slots = Int[]
        for a in args[2:end]
            au = Symbolics.unwrap(a)
            value = if au isa Number
                au
            else
                try
                    SymbolicUtils.unwrap_const(au)
                catch
                    return nothing
                end
            end
            value isa Integer || return nothing
            push!(slots, Int(value))
        end
        return isempty(slots) ? nothing : slots
    end
    return p isa SymbolicUtils.BasicSymbolic ? _param_slots(p) : nothing
end

function _kernel_parameter_sources(pmap)
    direct = Dict{Any, Any}()
    arrays = Dict{Symbol, Any}()
    named = Dict{Tuple{Symbol, Any}, Any}()
    for (key, value) in pmap
        ku = Symbolics.unwrap(key)
        direct[ku] = value
        name = _kernel_param_name(ku)
        name === nothing && continue
        slots = _kernel_param_slots(ku)
        if value isa AbstractArray
            arrays[name] = value
        else
            named[(name, slots === nothing ? nothing : Tuple(slots))] = value
        end
    end
    return direct, arrays, named
end

function _kernel_parameter_value(param, direct, arrays, named)
    haskey(direct, param) && return true, direct[param]
    name = _kernel_param_name(param)
    name === nothing && return false, nothing
    slots = _kernel_param_slots(param)
    slotkey = slots === nothing ? nothing : Tuple(slots)
    haskey(named, (name, slotkey)) && return true, named[(name, slotkey)]
    if slots !== nothing && haskey(arrays, name)
        return true, arrays[name][slots...]
    end
    return false, nothing
end

function _kernel_parameter_values(
        plan::KernelParameterPlan,
        pmap;
        current::Union{Nothing, Vector{ComplexF64}} = nothing,
        strict::Bool = true,
    )
    direct, arrays, named = _kernel_parameter_sources(pmap)
    values = current === nothing ? zeros(ComplexF64, length(plan.params)) : copy(current)
    missing = Any[]
    for (i, param) in enumerate(plan.params)
        found, value = _kernel_parameter_value(param, direct, arrays, named)
        if found
            try
                values[i] = ComplexF64(value)
            catch
                throw(
                    ArgumentError(
                        "kernel parameter $(param) must resolve to a numeric scalar; got " *
                            "$(typeof(value)).",
                    ),
                )
            end
        elseif current === nothing && strict
            push!(missing, param)
        end
    end
    isempty(missing) || throw(
        ArgumentError("missing values for kernel parameters: $(missing). Pass them in `ps`."),
    )
    return values
end

function _kernel_coefficient_values(plan::KernelParameterPlan, values::Vector{ComplexF64})
    pd = Dict{Any, Any}(plan.params[i] => values[i] for i in eachindex(plan.params))
    out = Vector{ComplexF64}(undef, length(plan.coeffs))
    for (i, coeff) in enumerate(plan.coeffs)
        if !(coeff isa Number)
            for var in Symbolics.get_variables(coeff)
                u = Symbolics.unwrap(var)
                if SymbolicUtils.issym(u) &&
                        Base.nameof(u) === :im &&
                        SymbolicUtils.symtype(u) === Number
                    pd[u] = im
                end
            end
        end
        substituted = Symbolics.substitute(coeff, pd)
        value = SymbolicUtils.unwrap_const(substituted)
        out[i] = try
            ComplexF64(value)
        catch err
            simplified = Symbolics.simplify(substituted; expand = true)
            simplified_value = SymbolicUtils.unwrap_const(simplified)
            try
                ComplexF64(simplified_value)
            catch
                throw(err)
            end
        end
    end
    return out
end

function KernelParameters(plan::KernelParameterPlan, pmap)
    values = _kernel_parameter_values(plan, pmap)
    return KernelParameters(values, _kernel_coefficient_values(plan, values))
end

# ---- direct problem construction -----------------------------------------------------

function _kernel_jacobian_flag(jac)
    jac === false && return false
    (jac === true || jac === :analytic) && return true
    throw(
        ArgumentError(
            "KernelBackend provides only its exact analytic Jacobian; `jac` must be " *
                "`false`, `true`, or `:analytic`, got $(repr(jac)).",
        ),
    )
end

function _build_kernel_rhs(eqs::MeanfieldEquations, ps, ::KernelBackend; jac = false)
    dojac = _kernel_jacobian_flag(jac)
    ir = _lower_moment_ir(eqs)
    kernel = MomentKernel(ir, ComplexF64)
    plan = KernelParameterPlan(ir)
    values = parameter_map(eqs, ps)
    p = KernelParameters(plan, values)
    jacobian = dojac ? MomentJacobianKernel(ir, ComplexF64) : nothing
    return KernelRHS(kernel, plan, jacobian), p
end

function _kernel_ode_function(rhs::KernelRHS)
    rhs.jacobian === nothing && return SciMLBase.ODEFunction{true}(rhs)
    return SciMLBase.ODEFunction{true}(
        rhs;
        jac = rhs.jacobian,
        jac_prototype = copy(rhs.jacobian.prototype),
    )
end

function _kernel_u0(eqs, u0::AbstractVector{<:Number})
    length(u0) == length(eqs.states) || throw(
        DimensionMismatch(
            "u0 length $(length(u0)) does not match number of states $(length(eqs.states))",
        ),
    )
    return ComplexF64.(u0)
end

function _kernel_u0(eqs, u0::AbstractDict)
    reg = _state_registry(eqs)
    return ComplexF64[
        haskey(u0, eqs.states[i]) ? ComplexF64(u0[eqs.states[i]]) :
            haskey(u0, reg.vars[i]) ? ComplexF64(u0[reg.vars[i]]) : zero(ComplexF64) for
            i in eachindex(eqs.states)
    ]
end

_kernel_u0(eqs, state) = initial_values(eqs, state)

"""
    SciMLBase.ODEProblem(eqs, u0, tspan, ps; backend = KernelBackend(), jac = false, kwargs...)

Construct a deterministic `ODEProblem` without converting `eqs` to a ModelingToolkit
`System`. The `backend` keyword is intentionally required while this specialized path is new.
Set `jac = true` (or `:analytic`) only for holomorphic/unfolded closures.
"""
function SciMLBase.ODEProblem(
        eqs::MeanfieldEquations,
        u0,
        tspan,
        ps;
        backend::KernelBackend,
        jac = false,
        kwargs...,
    )
    rhs, p = _build_kernel_rhs(eqs, ps, backend; jac)
    f = _kernel_ode_function(rhs)
    return SciMLBase.ODEProblem(f, _kernel_u0(eqs, u0), tspan, p; kwargs...)
end

# ---- parameter updates ---------------------------------------------------------------

function _kernel_rhs(prob::SciMLBase.ODEProblem)
    inner = try
        prob.f.f
    catch
        nothing
    end
    inner isa KernelRHS || throw(
        ArgumentError("update_parameters! requires an ODEProblem built with KernelBackend()."),
    )
    return inner
end

"""
    update_parameters!(prob, pairs)

Refresh the numeric parameter/coefficient payload of a direct kernel problem without
relowering equations or rebuilding the evaluator. Partial updates retain unspecified values.
"""
function update_parameters!(prob::SciMLBase.ODEProblem, pairs)
    rhs = _kernel_rhs(prob)
    p = prob.p
    p isa KernelParameters{ComplexF64} || throw(
        ArgumentError("kernel ODEProblem has an incompatible parameter payload $(typeof(p))."),
    )
    values = _kernel_parameter_values(rhs.plan, pairs; current = p.values, strict = false)
    coeffs = _kernel_coefficient_values(rhs.plan, values)
    copyto!(p.values, values)
    copyto!(p.coeffs, coeffs)
    return prob
end

# ---- direct solution lookup ----------------------------------------------------------

function _is_kernel_solution(sol)
    p = try
        sol.prob.p
    catch
        nothing
    end
    return p isa KernelParameters
end

function _kernel_state_match(avg, eqs::MeanfieldEquations)
    ctx = build_ctx(eqs)
    treatments = _treatments(eqs, ctx)
    ops = QAdd[(o = undo_average(s); o isa QAdd ? o : o * 1) for s in eqs.states]
    moments = MomentMap(ctx, treatments, ops, collect(Int32, 1:length(ops)))
    op = undo_average(avg)
    op isa QAdd || return nothing
    return match_moment(moments, op)
end

function get_solution(sol, avg::SymbolicUtils.BasicSymbolic, eqs::MeanfieldEquations)
    if !_is_kernel_solution(sol)
        return invoke(
            get_solution,
            Tuple{Any, SymbolicUtils.BasicSymbolic, AbstractMeanfieldEquations},
            sol,
            avg,
            eqs,
        )
    end

    r = _kernel_state_match(avg, eqs)
    r === nothing && throw(KeyError(avg))
    index, same = r
    i = Int(index)
    return same ? (τ -> _eval_at(sol, i, τ)) : (τ -> conj.(_eval_at(sol, i, τ)))
end
