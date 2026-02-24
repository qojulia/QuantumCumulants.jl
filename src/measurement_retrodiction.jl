function split_equations(eqin::NoiseEquations)::Tuple{MeanfieldEquations,MeanfieldEquations}
    determ = MeanfieldEquations(
        eqin.equations,
        eqin.operator_equations,
        eqin.states,
        eqin.operators,
        eqin.hamiltonian,
        eqin.jumps,
        eqin.jumps_dagger,
        eqin.rates,
        eqin.iv,
        eqin.varmap,
        eqin.order,
    )
    noise = MeanfieldEquations(
        eqin.noise_equations,
        eqin.operator_noise_equations,
        eqin.states,
        eqin.operators,
        eqin.hamiltonian,
        eqin.jumps,
        eqin.jumps_dagger,
        eqin.efficiencies, # 
        eqin.iv,
        eqin.varmap,
        eqin.order,
    )
    return determ, noise
end

function merge_backward_equations(
    determ::MeanfieldEquations,
    noise::MeanfieldEquations,
)::BackwardMeanfieldNoiseEquations
    return MeanfieldNoiseEquations(
        determ.equations,
        determ.operator_equations,
        noise.equations,
        noise.operator_equations,
        determ.states,
        determ.operators,
        determ.hamiltonian,
        determ.jumps,
        determ.jumps_dagger,
        determ.rates,
        noise.rates,
        determ.iv,
        determ.varmap,
        determ.order,
    )
end

# TODO: fix scale for MeanfieldNoiseEquations first
# function scale(he::BackwardMeanfieldNoiseEquations; kwargs...)
#     determ, noise = split_equations(he)
#     return merge_backward_equations(scale(determ), scale(noise))
# end

function meanfield_backward(
    a::Vector,
    H,
    J;
    Jdagger::Vector = adjoint.(J),
    rates = ones(Int, length(J)),
    efficiencies = zeros(Int, length(J)),
    multithread = false,
    simplify = true,
    order = nothing,
    mix_choice = maximum,
    iv = MTK.t_nounits,
)

    inds = vcat(get_indices(a), get_indices(H), get_indices(J))
    if isempty(inds)
        error("The function meanfield_backward() does not yet support indices.")
    end

    if rates isa Matrix
        error("Nondiagonal measurements are not supported")
    end
    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    rhs_noise = Vector{Any}(undef, length(a))
    imH = im*H

    JJd_JdJ_term =
        -sum(
            rates[i]*(average(J[i]*Jdagger[i]) - average(Jdagger[i]*J[i])) for
            i = 1:length(J)
        )

    function calculate_term!(i)
        rhs_ = commutator(-imH, a[i]) # backward: -H
        rhs_diss = _master_lindblad_backward(a[i], J, Jdagger, rates) # backward: recycling term J→J⁺ 
        rhs_trace = a[i]*JJd_JdJ_term # trace preserving term
        rhs_dY_dS = _dY_dS_extra_term(a[i], J, Jdagger, efficiencies .* rates) # extra term to be able to use Y(t)
        rhs[i] = rhs_ + rhs_diss + rhs_trace + rhs_dY_dS # sum

        rhs_noise[i] =
            _master_noise(a[i], adjoint(J), adjoint(Jdagger), efficiencies .* rates) # backward: J→J⁺
    end

    if multithread
        Threads.@threads for i = 1:length(a)
            calculate_term!(i)
        end
    else
        for i = 1:length(a)
            calculate_term!(i)
        end
    end

    # Average
    vs = map(average, a)
    rhs_avg = map(average, rhs)
    rhs_noise_avg = map(average, rhs_noise)

    if simplify
        rhs_avg = map(SymbolicUtils.simplify, rhs_avg)
        rhs_noise_avg = map(SymbolicUtils.simplify, rhs_noise_avg)
    end

    rhs = map(undo_average, rhs_avg)
    rhs_noise = map(undo_average, rhs_noise_avg) # this is not correct, see _master_noise (but also not used)

    if order !== nothing
        rhs_avg = [
            cumulant_expansion(r, order; simplify = simplify, mix_choice = mix_choice)
            for r ∈ rhs_avg
        ]
        rhs_noise_avg = [
            cumulant_expansion(r, order; simplify = simplify, mix_choice = mix_choice)
            for r ∈ rhs_noise_avg
        ]
    end

    eqs_avg = [Symbolics.Equation(l, r) for (l, r) in zip(vs, rhs_avg)]
    eqs = [Symbolics.Equation(l, r) for (l, r) in zip(a, rhs)]
    eqs_noise_avg = [Symbolics.Equation(l, r) for (l, r) in zip(vs, rhs_noise_avg)]
    eqs_noise = [Symbolics.Equation(l, r) for (l, r) in zip(a, rhs_noise)]
    varmap = make_varmap(vs, iv)

    me = BackwardMeanfieldNoiseEquations(
        eqs_avg,
        eqs,
        eqs_noise_avg,
        eqs_noise,
        vs,
        a,
        H,
        J,
        Jdagger,
        rates,
        efficiencies,
        iv,
        varmap,
        order,
    )
    return me
end

function _dY_dS_extra_term(a_, J, Jdagger, rates)
    args = Any[]
    for k = 1:length(J)
        if isequal(rates[k], 0)
            continue
        end
        if isa(rates[k], SymbolicUtils.Symbolic) ||
           isa(rates[k], Number) ||
           isa(rates[k], Function)
            # simplify on 
            c1 = rates[k]*average((J[k]*a_-average(J[k])*average(a_)))
            c2 = rates[k]*average((a_*Jdagger[k]-average(a_)*average(Jdagger[k])))
            SQA.push_or_append_nz_args!(args, -(c1+c2)*(average(Jdagger[k]+J[k])))
        elseif isa(rates[k], Matrix)
            error("Nondiagonal measurements are not supported")
        else
            error("Unknown rates type!")
        end
    end
    isempty(args) && return 0
    return sum(args)
end


function _master_lindblad_backward(a_, J, Jdagger, rates)
    args = Any[]
    for k = 1:length(J)
        if isa(rates[k], SymbolicUtils.Symbolic) ||
           isa(rates[k], Number) ||
           isa(rates[k], Function)
            c1 = -0.5*rates[k]*a_*Jdagger[k]*J[k]
            c2 = -0.5*rates[k]*Jdagger[k]*J[k]*a_
            c3 = rates[k]*J[k]*a_*Jdagger[k]
            SQA.push_or_append_nz_args!(args, c1)
            SQA.push_or_append_nz_args!(args, c2)
            SQA.push_or_append_nz_args!(args, c3)
        elseif isa(rates[k], Matrix)
            error("Nondiagonal measurements are not supported") # _master_lindblad_backward() is only used for retrodiction
        else
            error("Unknown rates type!")
        end
    end
    isempty(args) && return 0
    return QAdd(args)
end

function _append!(lhs::NoiseEquations, rhs::NoiseEquations)
    append!(lhs.noise_equations, rhs.noise_equations)
    append!(lhs.operator_noise_equations, rhs.operator_noise_equations)
    append!(lhs.equations, rhs.equations)
    append!(lhs.operator_equations, rhs.operator_equations)
    append!(lhs.states, rhs.states)
    append!(lhs.operators, rhs.operators)
    append!(lhs.varmap, rhs.varmap)
end

function MTK.complete!(
    de::BackwardMeanfieldNoiseEquations;
    order = de.order,
    multithread = false,
    filter_func = nothing,
    mix_choice = maximum,
    simplify = true,
    kwargs...,
)
    order = calculate_order(de, de.equations, order)
    order_noise = calculate_order(de, de.noise_equations, order)
    order = max(order, order_noise)
    missed = missing_variables(
        de,
        de.equations,
        order,
        multithread,
        filter_func,
        mix_choice,
        simplify,
    )
    missed_noise = missing_variables(
        de,
        de.noise_equations,
        order,
        multithread,
        filter_func,
        mix_choice,
        simplify,
    )
    missed = Set(vcat(missed, missed_noise))

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        he = meanfield_backward(
            ops_,
            de.hamiltonian,
            de.jumps;
            Jdagger = de.jumps_dagger,
            rates = de.rates,
            efficiencies = de.efficiencies,
            simplify = simplify,
            multithread = multithread,
            order = order,
            mix_choice = mix_choice,
            iv = de.iv,
            kwargs...,
        )
        _append!(de, he)
        missed = missing_variables(
            de,
            de.equations,
            order,
            multithread,
            filter_func,
            mix_choice,
            simplify,
        )
        missed_noise = missing_variables(
            de,
            de.noise_equations,
            order,
            multithread,
            filter_func,
            mix_choice,
            simplify,
        )
        missed = Set(vcat(missed, missed_noise))
    end
    # TODO: setting filtered expression to 0 is missing here!
    return de
end

function cumulant_expansion(
    de::BackwardMeanfieldNoiseEquations,
    order;
    multithread = false,
    mix_choice = maximum,
    kwargs...,
)
    determ, noise = split_equations(de)
    return merge_backward_equations(
        cumulant_expansion(determ, order; multithread, mix_choice, kwargs...),
        cumulant_expansion(noise, order; multithread, mix_choice, kwargs...),
    )
end

function MTK.complete(de::NoiseEquations; kwargs...)
    de_ = deepcopy(de)
    complete!(de_; kwargs...)
    return de_
end
