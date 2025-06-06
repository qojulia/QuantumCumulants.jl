# get_indices functions
get_indices(x::AvgSums) = get_indices(arguments(x))
get_indices(term::Average) = get_indices(arguments(term)[1])
function get_indices(term::QMul)
    args_nc = filter(x -> x isa IndexedOperator, term.args_nc)
    return unique(vcat(get_indices(args_nc),get_indices(term.arg_c)))
end
get_indices(a::IndexedOperator) = [a.ind]
get_indices(vec::AbstractVector) = unique(vcat(get_indices.(vec)...))
function get_indices(a::BasicSymbolic{DoubleIndexedVariable})
    meta = SymbolicUtils.metadata(a)[DoubleIndexedVariable]
    return unique([meta.ind1,meta.ind2])
end
get_indices(a::BasicSymbolic{IndexedVariable}) = [SymbolicUtils.metadata(a)[IndexedVariable].ind]
const Sums = Union{SingleSum,DoubleSum}
# get_indices(x::Sums) = unique(get_indices(arguments(x)))
get_indices(x::SingleSum) = get_indices(x.term)
get_indices(x::DoubleSum) = get_indices(x.innerSum.term)
get_indices(x::Number) = []
get_indices(term) = iscall(term) ? get_indices(arguments(term)) : []

#Usability functions:
Σ(a,b) = DoubleSum(a,b)  #Double-Sum here, because if variable a is not a single sum it will create a single sum anyway
Σ(a,b,c;kwargs...) = DoubleSum(a,b,c;kwargs...)
∑(args...; kwargs...) = Σ(args...; kwargs...)

IndexedOperator(x::IndexableOps,numb::Int64) = NumberedOperator(x,numb)
function IndexedOperator(x::QTerm, numb::Int64) # σ(1,1,2)
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    f([NumberedOperator(arg, numb) for arg in args]...)
end
IndexedVariable(x,numb::Int64) = SingleNumberedVariable(x,numb)
IndexedVariable(x,num1::Int64,num2::Int64;kwargs...) = DoubleNumberedVariable(x,num1,num2;kwargs...)
IndexedVariable(name::Symbol,ind1::Index,ind2::Index;kwargs...) = DoubleIndexedVariable(name,ind1,ind2;kwargs...)
DoubleIndexedVariable(x,num1::Int64,num2::Int64;kwargs...) = DoubleNumberedVariable(x,num1,num2;kwargs...)

#Numeric Conversion of NumberedOperators
function to_numeric(op::NumberedOperator,b::QuantumOpticsBase.CompositeBasis; ranges::Vector{Int64}=Int64[],kwargs...)
    if isempty(ranges)
        error("When calling to_numeric for indexed Operators, specification of the \"ranges\" keyword is needed! This keyword requires a vector of Integers, which specify the maximum range of the index for each hilbertspace.")
    end
    h = hilbert(op)
    if h isa ProductSpace
        if length(h.spaces) != length(ranges)
            error("Unequal length of hilbertspaces and ranges!")
        end
    else
        if length(ranges) != 1
            error("Wrong number of entries in ranges!")
        end
    end
    start = 0
    if h !== nothing #this is fine here since there are assertions above
        aon_ = acts_on(op)
        for i = 1:(aon_ - 1)
            start = start + ranges[i]
        end
    end
    aon = 0
    if start == 0
        aon = getNumber(op)[1] - 1
    else
        aon = op.numb + start
    end
    op_ = _to_numeric(op.op,b.bases[aon];kwargs...)
    return QuantumOpticsBase.LazyTensor(b,[aon],(op_,))
end
#function that returns the conjugate of an average, but also preserving the correct ordering
function _inconj(v::Average)
    f = operation(v)
    if f == conj
        return _inconj(arguments(v)[1])
    end
    arg = v.arguments[1]
    adj_arg = inadjoint(arg)
    return _average(adj_arg)
end
function _inconj(v::T) where T <: SymbolicUtils.BasicSymbolic
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(_inconj, SymbolicUtils.arguments(v))
        return SymbolicUtils.maketerm(T, f, args, TermInterface.metadata(v))
    else
        return conj(v)
    end
end
_inconj(x::Number) = conj(x)

function inadjoint(q::QMul)
    qad = adjoint(q)
    inorder!(qad)
    return qad
end
inadjoint(op::QNumber) = adjoint(op)
inadjoint(s::SymbolicUtils.BasicSymbolic{<:Number}) = _conj(s)
inadjoint(x) = adjoint(x)

function inorder!(v::Average)
    f = operation(v)
    if f == conj
        return conj(inorder!(arguments(v)[1]))
    end
    return average(inorder!(arguments(v)[1]))
end
function inorder!(q::QMul)
    sort!(q.args_nc, by=get_numbers)
    sort!(q.args_nc, by=getIndName)
    sort!(q.args_nc, by=acts_on)
    return merge_commutators(q.arg_c,q.args_nc)
end
function inorder!(v::T) where T <: SymbolicUtils.BasicSymbolic
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(inorder!, SymbolicUtils.arguments(v))
        return SymbolicUtils.maketerm(T, f, args, TermInterface.metadata(v))
    end
    return v
end
inorder!(x) = x

get_numbers(term::Average) = get_numbers(arguments(term)[1])
get_numbers(term::QMul) = unique(vcat(get_numbers.(term.args_nc)...))
get_numbers(x::NumberedOperator) = [x.numb]
get_numbers(x::AbstractVector) = unique(vcat(get_numbers.(x)...))
get_numbers(x) = []


"""
    subst_reds(de::AbstractMeanfieldEquations)

Function that substitutes possible redundant conjugate averages inside the given Equations with their corresponding
average given as the conjugate of one of the left-hand-side (of the equations) averages.

# Optional Arguments
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.
"""
subst_reds(de::AbstractMeanfieldEquations;scaling=false,kwargs...) = scaling ? subst_reds_scale(de;kwargs...) : subst_reds_eval(de;kwargs...)

function subst_reds_eval(me::AbstractMeanfieldEquations;kwargs...)
    states = deepcopy(me.states)
    to_sub = _inconj.(states)
    filter!(x->x ∉ me.states, to_sub)
    to_insert = conj.(_inconj.(to_sub))
    subs = Dict(to_sub .=> to_insert)
    eqs = [substitute(eq,subs) for eq in me.equations]
    return IndexedMeanfieldEquations(eqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end

function subst_reds_scale(me::AbstractMeanfieldEquations;kwargs...)

    states = me.states
    missed = find_missing(me)
    inorder!.(missed)

    to_sub = missed
    filter!(x->x ∉ states,to_sub)

    to_insert = Vector{Any}(nothing,length(to_sub))

    counter = 1
    while counter <= length(to_sub)
        elem = to_sub[counter]
        ind_ = findfirst(x -> isscaleequal(elem,x;kwargs...),states)
        if !=(ind_,nothing)
            to_insert[counter] = states[ind_]
            counter = counter + 1
        else
            ind_ = findfirst(x -> isscaleequal(inorder!(_inconj(elem)),x;kwargs...),states)
            if !=(ind_,nothing)
                to_insert[counter] = conj(states[ind_])
                counter = counter + 1
            else
                deleteat!(to_insert,counter) # these deletes are for consistancy only -> it is possible that not all terms are fully evaluated
                deleteat!(to_sub,counter)   # yet in the system -> leftovers in the find_missing
            end
        end
    end

    subs = Dict(to_sub .=> to_insert)
    eqs = [substitute(eq,subs) for eq in me.equations]

    return IndexedMeanfieldEquations(eqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end
function subst_reds_scale(term::SymbolicUtils.BasicSymbolic;kwargs...)
    avrgs = unique(getAvrgs(term))
    D = Dict{Average,Average}()
    for i = 1:length(avrgs)
        ind_ = findfirst(x -> isscaleequal(avrgs[i],x;kwargs...) && !isequal(avrgs[i],x),avrgs)
        if !=(ind_,nothing)
            push!(D,avrgs[i] => avrgs[ind_])
            avrgs[i] = nothing
        end
    end
    return inorder!(substitute(term,D;kwargs...))
end

#function that checks if 2 averages are the same, if they would get scaled
# or when they are within a scaled system
# meaning, that for example <σ²²₁ σ²¹₂> == <σ²¹₁ σ²²₂> when both σ are acting on the same hilbertspace
isscaleequal(avr1::Average,avr2::Average;kwargs...) = isscaleequal(arguments(avr1)[1],arguments(avr2)[1];kwargs...)
function isscaleequal(qmul1::QMul,qmul2::QMul;kwargs...)
    isequal(qmul1,qmul2) && return true
    isequal(length(qmul1.args_nc), length(qmul2.args_nc)) || return false
    isequal(acts_on(qmul1),acts_on(qmul2)) || return false
    return has_same(get_ops_of_aons(qmul1;kwargs...),get_ops_of_aons(qmul2;kwargs...))
end
function isscaleequal(a::NumberedOperator,b::NumberedOperator;h=nothing,kwargs...)
    isequal(a,b) && return true
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        acts_on(a) == acts_on(b) || return false
        (acts_on(a) in h) && return isequal(a.op,b.op)
        return isequal(a,b)
    else
        return isequal(a.op,b.op)
    end
end
function isscaleequal(a::IndexedOperator,b::IndexedOperator;h=nothing,kwargs...)
    isequal(a,b) && return true
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        acts_on(a) == acts_on(b) || return false
        (acts_on(a) in h) && return isequal(a.op,b.op)
        return isequal(a,b)
    else
        return isequal(a.op,b.op)
    end
end
function isscaleequal(t1,t2;kwargs...)
    if SymbolicUtils.iscall(t1) && SymbolicUtils.iscall(t2)
        args1 = arguments(t1)
        args2 = arguments(t2)
        isequal(operation(t1),operation(t2)) || return false
        length(args1) != length(args2) && return false
        for i = 1:length(args1)
            isscaleequal(args1[i],args2[i];kwargs...) || return false
        end
        return true
    end
    return isequal(t1,t2)
end
function has_same(vec1::Vector,vec2::Vector)
    length(vec1) != length(vec2) && return false
    return isequal(counter.(vec1),counter.(vec2))
end
#function that returns a vector of vector of ops depending on their aons
# so all ops, that act on the same  are returned in the same subvector
function get_ops_of_aons(qmul::QMul; h=nothing, kwargs...)
    aons = acts_on(qmul)
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        filter!(x -> x in h,aons)
    end
    dic = Dict{Int,Any}(aon => Any[] for aon in aons)
    noAonVec = []
    for op in qmul.args_nc
        if acts_on(op) in aons
            if (op isa NumberedOperator || op isa IndexedOperator) && op.op isa Transition
                push!(dic[acts_on(op)],op.op)
            else
                push!(dic[acts_on(op)],op)
            end
        else
            push!(noAonVec,op)
        end
    end
    allVec = []
    for i in keys(dic)
        push!(allVec,dic[i])
    end
    push!(allVec,noAonVec)
    return allVec
end

# function, that creates a dictionary within each key is an element of it and the value the count of the elemt,
# used to calculate, if two arrays are containing the exact same elements, without order
# copied from https://discourse.julialang.org/t/compare-array-of-string-with-no-regard-to-the-order/62322
function counter(it)
    y = Dict{eltype(it), Int}()
    for i in it
        y[i] = get(y, i, 0) + 1
    end
    return y
end
function isNotIn(avrg::Average,states::Vector,scaling;kwargs...)
    _avg = inorder!(avrg)
    if scaling
        for state in states
            if isscaleequal(avrg,state;kwargs...) || isscaleequal(_inconj(avrg),state;kwargs...)
                return false
            end
        end
        return true
    else
        for state in states # i tried to do something like _avg ∉ states, this somehow interferes with the order of operators in the average
            isequal(_avg,state) && return false
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
function value_map(ps::Vector,p0::Vector;limits=nothing,kwargs...)
    length(ps) != length(p0) && error("Vectors given have non-equal length!")

    if !=(limits,nothing) && limits isa Pair
        mapping_ = Dict{BasicSymbolic,Int64}(first(limits)=>last(limits))
        limits = mapping_
    end
    if limits === nothing
        limits = Dict{BasicSymbolic,Int64}()
    end

    dict = Dict{SymbolicUtils.BasicSymbolic,ComplexF64}()
    for i=1:length(ps)
        dicVal = nothing
        if ps[i] isa BasicSymbolic{IndexedVariable}
            if p0[i] isa Vector || p0[i] isa Number
                dicVal = create_value_map(ps[i],p0[i];limits)
            else
                error("cannot resolve entry at $i-th position in values-vector")
            end
        elseif ps[i] isa BasicSymbolic{DoubleIndexedVariable}
            if p0[i] isa Matrix || p0[i] isa Number
                dicVal = create_value_map(ps[i],p0[i];limits)
            end
        else
            push!(dict,ps[i]=>p0[i])
            continue
        end
        dict = merge(dict,dicVal)
    end
    return collect(dict)
end
