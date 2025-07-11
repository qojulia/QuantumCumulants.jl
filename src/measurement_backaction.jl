"""
    MeanfieldNoiseEquations

Mean field equations including a separate set of equations describing the
noise generated by measurement backactions.
"""
struct MeanfieldNoiseEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation} # useless but needed to create eqs
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    efficiencies::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

function _master_noise(a_, J, Jdagger, rates)
    args = Any[]
    for k = 1:length(J)
        if isequal(rates[k], 0)
            continue
        end
        if isa(rates[k], SymbolicUtils.Symbolic) ||
           isa(rates[k], Number) ||
           isa(rates[k], Function)
            c1 = sqrt(rates[k])*(Jdagger[k]*a_-average(Jdagger[k])*average(a_))
            c2 = sqrt(rates[k])*(a_*J[k]-average(a_)*average(J[k]))
            SQA.push_or_append_nz_args!(args, c1)
            SQA.push_or_append_nz_args!(args, c2)
        elseif isa(rates[k], Matrix)
            error("Nondiagonal measurements are not supported")
        else
            error("Unknown rates type!")
        end
    end
    isempty(args) && return 0
    return QAdd(args)
end

function split_equations(
    eqin::MeanfieldNoiseEquations,
)::Tuple{MeanfieldEquations,MeanfieldEquations}
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
        eqin.efficiencies,
        eqin.iv,
        eqin.varmap,
        eqin.order,
    )
    return determ, noise
end

function merge_equations(
    determ::MeanfieldEquations,
    noise::MeanfieldEquations,
)::MeanfieldNoiseEquations
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

function scale(he::MeanfieldNoiseEquations; kwargs...)
    determ, noise = split_equations(he)
    return merge_equations(scale(determ), scale(noise))
end

function _meanfield_backaction(
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

    if rates isa Matrix
        J = [J];
        Jdagger = [Jdagger];
        rates = [rates]
    end
    J_, Jdagger_, rates_ = _expand_clusters(J, Jdagger, rates)
    J_, Jdagger_, efficiencies_ = _expand_clusters(J, Jdagger, efficiencies)
    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    rhs_noise = Vector{Any}(undef, length(a))
    imH = im*H

    function calculate_term(i)
        rhs_ = commutator(imH, a[i])
        rhs_diss = _master_lindblad(a[i], J_, Jdagger_, rates_)
        rhs_noise[i] = _master_noise(a[i], J_, Jdagger_, efficiencies_ .* rates)
        rhs[i] = rhs_ + rhs_diss
    end

    if multithread
        Threads.@threads for i = 1:length(a)
            calculate_term(i)
        end
    else
        for i = 1:length(a)
            calculate_term(i)
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
    rhs_noise = map(undo_average, rhs_noise_avg)

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

    me = MeanfieldNoiseEquations(
        eqs_avg,
        eqs,
        eqs_noise_avg,
        eqs_noise,
        vs,
        a,
        H,
        J_,
        Jdagger_,
        rates_,
        efficiencies_,
        iv,
        varmap,
        order,
    )
    if has_cluster(H)
        return scale(me; simplify = simplify, order = order, mix_choice = mix_choice)
    else
        return me
    end
end


function calculate_order(de::AbstractMeanfieldEquations, eqns, order)
    vs = de.states
    order_lhs = maximum(get_order.(vs))
    order_rhs = 0
    for i = 1:length(eqns)
        k = get_order(eqns[i].rhs)
        k > order_rhs && (order_rhs = k)
    end
    if order === nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end

    maximum(order_) >= order_lhs || error(
        "Cannot form cumulant expansion of derivative; you may want to use a higher order!",
    )
    return order_
end


function missing_variables(
    de::AbstractMeanfieldEquations,
    eqns,
    order = de.order,
    multithread = false,
    filter_func = nothing,
    mix_choice = maximum,
    simplify = true,
)
    vs = de.states
    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(eqns, vhash, vs′hash; get_adjoints = false)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    return missed
end


function filter_set_zero!(
    de::AbstractMeanfieldEquations,
    order = de.order,
    multithread = false,
    filter_func = nothing,
    mix_choice = maximum,
    simplify = true,
)
    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints = false)
        filter!(!filter_func, missed)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i = 1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs)
            de.states[i] = de.equations[i].lhs
        end
    end
end

function _append!(lhs::MeanfieldNoiseEquations, rhs::MeanfieldNoiseEquations)
    append!(lhs.noise_equations, rhs.noise_equations)
    append!(lhs.operator_noise_equations, rhs.operator_noise_equations)
    append!(lhs.equations, rhs.equations)
    append!(lhs.operator_equations, rhs.operator_equations)
    append!(lhs.states, rhs.states)
    append!(lhs.operators, rhs.operators)
    append!(lhs.varmap, rhs.varmap)
end

function MTK.complete!(
    de::MeanfieldNoiseEquations;
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
        he = _meanfield_backaction(
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
    return de
end

function cumulant_expansion(
    de::MeanfieldNoiseEquations,
    order;
    multithread = false,
    mix_choice = maximum,
    kwargs...,
)
    determ, noise = split_equations(de)
    return merge_equations(
        cumulant_expansion(determ, order; multithread, mix_choice, kwargs...),
        cumulant_expansion(noise, order; multithread, mix_choice, kwargs...),
    )
end


function MTK.complete(de::MeanfieldNoiseEquations; kwargs...)
    de_ = deepcopy(de)
    complete!(de_; kwargs...)
    return de_
end
