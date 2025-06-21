"""
    subst_reds(de::AbstractMeanfieldEquations)

Function that substitutes possible redundant conjugate averages inside the given Equations with their corresponding
average given as the conjugate of one of the left-hand-side (of the equations) averages.

# Optional Arguments
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.
"""
subst_reds(de::AbstractMeanfieldEquations; scaling = false, kwargs...) =
    scaling ? subst_reds_scale(de; kwargs...) : subst_reds_eval(de; kwargs...)

function subst_reds_eval(me::AbstractMeanfieldEquations; kwargs...)
    states = deepcopy(me.states)
    to_sub = _inconj.(states)
    filter!(x->x ∉ me.states, to_sub)
    to_insert = conj.(_inconj.(to_sub))
    subs = Dict(to_sub .=> to_insert)
    eqs = [substitute(eq, subs) for eq in me.equations]
    return IndexedMeanfieldEquations(
        eqs,
        me.operator_equations,
        me.states,
        me.operators,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.iv,
        me.varmap,
        me.order,
    )
end

function subst_reds_scale(me::AbstractMeanfieldEquations; kwargs...)

    states = me.states
    missed = find_missing(me)
    inorder!.(missed)

    to_sub = missed
    filter!(x->x ∉ states, to_sub)

    to_insert = Vector{Any}(nothing, length(to_sub))

    counter = 1
    while counter <= length(to_sub)
        elem = to_sub[counter]
        ind_ = findfirst(x -> isscaleequal(elem, x; kwargs...), states)
        if !=(ind_, nothing)
            to_insert[counter] = states[ind_]
            counter = counter + 1
        else
            ind_ =
                findfirst(x -> isscaleequal(inorder!(_inconj(elem)), x; kwargs...), states)
            if !=(ind_, nothing)
                to_insert[counter] = conj(states[ind_])
                counter = counter + 1
            else
                deleteat!(to_insert, counter) # these deletes are for consistancy only -> it is possible that not all terms are fully evaluated
                deleteat!(to_sub, counter)   # yet in the system -> leftovers in the find_missing
            end
        end
    end

    subs = Dict(to_sub .=> to_insert)
    eqs = [substitute(eq, subs) for eq in me.equations]

    return IndexedMeanfieldEquations(
        eqs,
        me.operator_equations,
        me.states,
        me.operators,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.iv,
        me.varmap,
        me.order,
    )
end
function subst_reds_scale(term::SymbolicUtils.BasicSymbolic; kwargs...)
    avrgs = unique(getAvrgs(term))
    D = Dict{Average,Average}()
    for i = 1:length(avrgs)
        ind_ = findfirst(
            x -> isscaleequal(avrgs[i], x; kwargs...) && !isequal(avrgs[i], x),
            avrgs,
        )
        if !=(ind_, nothing)
            push!(D, avrgs[i] => avrgs[ind_])
            avrgs[i] = nothing
        end
    end
    return inorder!(substitute(term, D; kwargs...))
end

#function that checks if 2 averages are the same, if they would get scaled
# or when they are within a scaled system
# meaning, that for example <σ²²₁ σ²¹₂> == <σ²¹₁ σ²²₂> when both σ are acting on the same hilbertspace
isscaleequal(avr1::Average, avr2::Average; kwargs...) =
    isscaleequal(arguments(avr1)[1], arguments(avr2)[1]; kwargs...)
function isscaleequal(qmul1::QMul, qmul2::QMul; kwargs...)
    isequal(qmul1, qmul2) && return true
    isequal(length(qmul1.args_nc), length(qmul2.args_nc)) || return false
    isequal(acts_on(qmul1), acts_on(qmul2)) || return false
    return has_same(get_ops_of_aons(qmul1; kwargs...), get_ops_of_aons(qmul2; kwargs...))
end
function isscaleequal(a::NumberedOperator, b::NumberedOperator; h = nothing, kwargs...)
    isequal(a, b) && return true
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        acts_on(a) == acts_on(b) || return false
        (acts_on(a) in h) && return isequal(a.op, b.op)
        return isequal(a, b)
    else
        return isequal(a.op, b.op)
    end
end
function isscaleequal(a::IndexedOperator, b::IndexedOperator; h = nothing, kwargs...)
    isequal(a, b) && return true
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        acts_on(a) == acts_on(b) || return false
        (acts_on(a) in h) && return isequal(a.op, b.op)
        return isequal(a, b)
    else
        return isequal(a.op, b.op)
    end
end
function isscaleequal(t1, t2; kwargs...)
    if SymbolicUtils.iscall(t1) && SymbolicUtils.iscall(t2)
        args1 = arguments(t1)
        args2 = arguments(t2)
        isequal(operation(t1), operation(t2)) || return false
        length(args1) != length(args2) && return false
        for i = 1:length(args1)
            isscaleequal(args1[i], args2[i]; kwargs...) || return false
        end
        return true
    end
    return isequal(t1, t2)
end
function has_same(vec1::Vector, vec2::Vector)
    length(vec1) != length(vec2) && return false
    return isequal(counter.(vec1), counter.(vec2))
end


# function, that creates a dictionary within each key is an element of it and the value the count of the elemt,
# used to calculate, if two arrays are containing the exact same elements, without order
# copied from https://discourse.julialang.org/t/compare-array-of-string-with-no-regard-to-the-order/62322
function counter(it)
    y = Dict{eltype(it),Int}()
    for i in it
        y[i] = get(y, i, 0) + 1
    end
    return y
end
function isNotIn(avrg::Average, states::Vector, scaling; kwargs...)
    _avg = inorder!(avrg)
    if scaling
        for state in states
            if isscaleequal(avrg, state; kwargs...) ||
               isscaleequal(_inconj(avrg), state; kwargs...)
                return false
            end
        end
        return true
    else
        for state in states # i tried to do something like _avg ∉ states, this somehow interferes with the order of operators in the average
            isequal(_avg, state) && return false
        end
        return true
    end
end

"""
    value_map(ps::Vector,p0::Vector)

A Function to create parameter values for indexed Variables more convenient.

# Arguments
*`ps::Vector`: A vector of parameters, that have no value assigned to them.
*`p0::Vector`: A vector for numeric values, that should get assigned to the corresponding
    entry in the `ps` vector. For Single-Indexed Variables the entry in the vector can also be again
    a Vector, that has an amount of entries as the index of the variables has range. For Double-Indexed
    Variables, this can also be a Matrix of a dimension, that corresponds to the ranges of the indices
    of the given variable.

"""
function value_map(ps::Vector, p0::Vector; limits = nothing, kwargs...)
    length(ps) != length(p0) && error("Vectors given have non-equal length!")

    if !=(limits, nothing) && limits isa Pair
        mapping_ = Dict{BasicSymbolic,Int64}(first(limits)=>last(limits))
        limits = mapping_
    end
    if limits === nothing
        limits = Dict{BasicSymbolic,Int64}()
    end

    dict = Dict{SymbolicUtils.BasicSymbolic,ComplexF64}()
    for i = 1:length(ps)
        dicVal = nothing
        if ps[i] isa BasicSymbolic{IndexedVariable}
            if p0[i] isa Vector || p0[i] isa Number
                dicVal = create_value_map(ps[i], p0[i]; limits)
            else
                error("cannot resolve entry at $i-th position in values-vector")
            end
        elseif ps[i] isa BasicSymbolic{DoubleIndexedVariable}
            if p0[i] isa Matrix || p0[i] isa Number
                dicVal = create_value_map(ps[i], p0[i]; limits)
            end
        else
            push!(dict, ps[i]=>p0[i])
            continue
        end
        dict = merge(dict, dicVal)
    end
    return collect(dict)
end


#function that returns a vector of vector of ops depending on their aons
# so all ops, that act on the same  are returned in the same subvector
function get_ops_of_aons(qmul::QMul; h = nothing, kwargs...)
    aons = acts_on(qmul)
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        filter!(x -> x in h, aons)
    end
    dic = Dict{Int,Any}(aon => Any[] for aon in aons)
    noAonVec = []
    for op in qmul.args_nc
        if acts_on(op) in aons
            if (op isa NumberedOperator || op isa IndexedOperator) && op.op isa Transition
                push!(dic[acts_on(op)], op.op)
            else
                push!(dic[acts_on(op)], op)
            end
        else
            push!(noAonVec, op)
        end
    end
    allVec = []
    for i in keys(dic)
        push!(allVec, dic[i])
    end
    push!(allVec, noAonVec)
    return allVec
end
