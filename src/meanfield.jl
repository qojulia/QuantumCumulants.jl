"""
    meanfield(ops, H, J=QField[]; Jdagger=adjoint.(J), rates=ones(length(J)),
              efficiencies=nothing, direction=Forward(),
              order=nothing, simplify=true,
              mix_choice=maximum, iv=Symbolics.variable(:t))

Compute equations of motion for the averages of `ops` under Hamiltonian `H`
and collapse operators `J` (with rates `rates`). Returns a `MeanFieldEquations`
unless `efficiencies` is given (then a `NoiseMeanFieldEquations`).
"""
_make_iv() = first(MTK.@independent_variables t)

function meanfield(
    ops::AbstractVector,
    H::QField,
    J::AbstractVector = QField[];
    Jdagger::AbstractVector = adjoint.(J),
    rates::AbstractVector = ones(length(J)),
    efficiencies = nothing,
    direction::EvolutionDirection = Forward(),
    order = nothing,
    simplify::Bool = true,
    mix_choice::Function = maximum,
    iv::Symbolics.Num = _make_iv(),
)
    Jn, Jdn = _normalize_jumps(J, Jdagger)
    rn = _normalize_rates(rates, length(Jn))
    if efficiencies === nothing
        direction isa Forward || throw(ArgumentError(
            "`direction=Backward()` requires `efficiencies=...`"))
        return _meanfield_forward(ops, H, Jn, Jdn, rn, order,
                                  simplify, mix_choice, iv)
    else
        en = _normalize_rates(efficiencies, length(Jn))
        return _meanfield_noise(direction, ops, H, Jn, Jdn, rn, en,
                                order, simplify, mix_choice, iv)
    end
end

function _normalize_jumps(J, Jdagger)
    if isempty(J)
        return QField[], QField[]
    end
    return collect(J), collect(Jdagger)
end

function _normalize_rates(rates, n::Int)
    if isempty(rates) && n == 0
        return Symbolics.Num[]
    end
    return collect(rates)
end

meanfield(op::QField, H::QField, args...; kw...) = meanfield([op], H, args...; kw...)

function _meanfield_forward(ops, H, J, Jdagger, rates, order,
                            simplify, mix_choice, iv)
    imH = im * H
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for (i, op) in enumerate(ops_qa)
        rhs = commutator(imH, op) + _lindblad_rhs(op, J, Jdagger, rates)
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    states  = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
                _normalize_order(order, (; hamiltonian = H))
    avg_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i])
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; simplify=false, mix_choice)
        end
        simplify && (rhs = SymbolicUtils.simplify(rhs))
        avg_eqs[i] = states[i] ~ rhs
    end
    return MeanFieldEquations(avg_eqs, operator_eqs, states, ops_qa,
                              H, collect(J), collect(Jdagger), collect(rates),
                              iv, order_vec)
end

function _lindblad_rhs(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        acc += (rates[k]/2) * (Jdagger[k] * commutator(op, J[k]) +
                                commutator(Jdagger[k], op) * J[k])
    end
    return acc
end

function _meanfield_noise(::Forward, ops, H, J, Jdagger, rates,
                          efficiencies, order, simplify, mix_choice, iv)
    imH = im * H
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for (i, op) in enumerate(ops_qa)
        rhs = commutator(imH, op) + _lindblad_rhs(op, J, Jdagger, rates)
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
                _normalize_order(order, (; hamiltonian = H))
    avg_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i])
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; simplify=false, mix_choice)
        end
        simplify && (rhs = SymbolicUtils.simplify(rhs))
        avg_eqs[i] = states[i] ~ rhs
    end
    op_noise, avg_noise = _build_noise_equations_forward(
        ops_qa, J, Jdagger, rates, efficiencies, simplify)
    if order_vec !== nothing
        avg_noise = [eq.lhs ~ cumulant_expansion(eq.rhs, order_vec;
                                                  simplify=false, mix_choice)
                     for eq in avg_noise]
    end
    if simplify
        avg_noise = [eq.lhs ~ SymbolicUtils.simplify(eq.rhs) for eq in avg_noise]
    end
    return NoiseMeanFieldEquations(avg_eqs, avg_noise, operator_eqs, op_noise,
                                    states, ops_qa, H,
                                    collect(J), collect(Jdagger),
                                    collect(rates), collect(efficiencies),
                                    iv, order_vec, Forward())
end

function _meanfield_noise(::Backward, ops, H, J, Jdagger, rates,
                          efficiencies, order, simplify, mix_choice, iv)
    imH = im * H
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for (i, op) in enumerate(ops_qa)
        rhs = commutator(imH, op) + _lindblad_rhs(op, J, Jdagger, rates)
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
                _normalize_order(order, (; hamiltonian = H))
    avg_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i])
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; simplify=false, mix_choice)
        end
        simplify && (rhs = SymbolicUtils.simplify(rhs))
        avg_eqs[i] = states[i] ~ rhs
    end
    op_noise, avg_noise = _build_noise_equations_backward(
        ops_qa, J, Jdagger, rates, efficiencies, simplify)
    if order_vec !== nothing
        avg_noise = [eq.lhs ~ cumulant_expansion(eq.rhs, order_vec;
                                                  simplify=false, mix_choice)
                     for eq in avg_noise]
    end
    if simplify
        avg_noise = [eq.lhs ~ SymbolicUtils.simplify(eq.rhs) for eq in avg_noise]
    end
    return NoiseMeanFieldEquations(avg_eqs, avg_noise, operator_eqs, op_noise,
                                    states, ops_qa, H,
                                    collect(J), collect(Jdagger),
                                    collect(rates), collect(efficiencies),
                                    iv, order_vec, Backward())
end
