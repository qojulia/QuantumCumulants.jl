# get_indices functions
get_indices(x::AvgSums) = get_indices(arguments(x))
get_indices(term::Average) = get_indices(arguments(term)[1])
function get_indices(term::QMul)
    args_nc = filter(x -> x isa IndexedOperator, term.args_nc)
    return unique(vcat(get_indices(args_nc),get_indices(term.arg_c)))
end
get_indices(a::IndexedOperator) = [a.ind]
get_indices(vec::Vector) = unique(vcat(get_indices.(vec)...))
get_indices(a::SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}) = unique([a.metadata.ind1,a.metadata.ind2])
get_indices(a::SymbolicUtils.Sym{Parameter,IndexedVariable}) = [a.metadata.ind]
const Sums = Union{SingleSum,IndexedDoubleSum}
get_indices(x::Sums) = unique(get_indices(arguments(x)))
get_indices(x::Number) = []
get_indices(term) = istree(term) ? get_indices(arguments(term)) : []

#Usability functions:
Σ(a,b) = IndexedDoubleSum(a,b)  #Double-Sum here, because if variable a is not a single sum it will create a single sum anyway
Σ(a,b,c;kwargs...) = IndexedDoubleSum(a,b,c;kwargs...)
∑(args...; kwargs...) = Σ(args...; kwargs...)

IndexedOperator(x::IndexableOps,numb::Int64) = NumberedOperator(x,numb)
IndexedVariable(x,numb::Int64) = SingleNumberedVariable(x,numb)
IndexedVariable(x,num1::Int64,num2::Int64;kwargs...) = DoubleNumberedVariable(x,num1,num2;kwargs...)
IndexedVariable(name::Symbol,ind1::Index,ind2::Index;kwargs...) = DoubleIndexedVariable(name,ind1,ind2;kwargs...)

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
    if !=(h,nothing) #this is fine here since there are assertions above
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
    return QuantumOpticsBase.embed(b,aon,op_)
end
#function that returns the conjugate of an average, but also preserving the correct ordering
function _inconj(v::Average)
    arg = v.arguments[1]
    adj_arg = inadjoint(arg)
    return _average(adj_arg)
end
function inadjoint(q::QMul)
    qad = adjoint(q)
    inorder!(qad)
    return qad
end
inadjoint(op::QNumber) = adjoint(op)
inadjoint(s::SymbolicUtils.Symbolic{<:Number}) = _conj(s)
inadjoint(x) = adjoint(x)

inorder!(v::Average) = average(inorder!(arguments(v)[1]))
function inorder!(q::QMul)
    sort!(q.args_nc, by=_get_number)
    sort!(q.args_nc, by=getIndName)
    sort!(q.args_nc, by=acts_on)
    return q
end
inorder!(x) = x

get_numbers(term::Average) = get_numbers(arguments(term)[1])
get_numbers(term::QMul) = unique(vcat(_get_number.(term.args_nc)))
_get_number(x::NumberedOperator) = x.numb
_get_number(x) = 0


"""
    subst_reds(de::MeanfieldEquations)

Function that substitutes possible redundant conjugate averages inside the given Equations with their corresponding
average given as the conjugate of one of the left-hand-side (of the equations) averages.

# Optional Arguments
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.
"""
function subst_reds(me::AbstractMeanfieldEquations;scaling::Bool=false,kwargs...)
    to_sub = find_missing(me)
    to_insert = Vector{Any}(nothing,length(to_sub))
    to_sub = inorder!.(to_sub)
    filter!(x->!(x in me.states), to_sub) #this one might be redundant
    if scaling
        states = deepcopy(me.states)
        counter = 1
        while counter <= length(to_sub)
            elem = to_sub[counter]
            ind_ = findfirst(x -> isscaleequal(elem,x;kwargs...),states)
            if !=(ind_,nothing)
                to_insert[counter] = states[ind_]
                counter = counter + 1
                continue
            end
            ind_ = findfirst(x -> isscaleequal(_inconj(elem),x;kwargs...),states)
            if !=(ind_,nothing)
                to_insert[counter] = conj(states[ind_])
                counter = counter + 1
                continue
            end
            deleteat!(to_insert,counter)
            deleteat!(to_sub,counter)
        end
    else  
        to_insert = conj(_inconj.(to_sub))
    end
    filter!(x -> !=(x,nothing),to_insert)
    subs = Dict(to_sub .=> to_insert)
    eqs = [substitute(eq,subs) for eq in me.equations]
    return IndexedMeanfieldEquations(eqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end

#function that checks if 2 averages are the same, if they would get scaled
# or within a scaled system
# meaning, that for example <σ²²₁ σ²¹₂> == <σ²¹₁ σ²²₂>   
isscaleequal(avr1::Average,avr2::Average;kwargs...) = isscaleequal(arguments(avr1)[1],arguments(avr2)[1];kwargs...)
function isscaleequal(qmul1::QMul,qmul2::QMul;kwargs...)
    isequal(qmul1,qmul2) && return true
    isequal(length(qmul1.args_nc), length(qmul2.args_nc)) || return false

    return has_same(vcat(get_ops_of_inds(qmul1;kwargs...),get_ops_of_numbs(qmul1;kwargs...)),vcat(get_ops_of_inds(qmul2;kwargs...),get_ops_of_numbs(qmul2;kwargs...)))
end
isscaleequal(a::NumberedOperator,b::NumberedOperator;h=nothing,kwargs...) = isequal(a.op,b.op) 
isscaleequal(a::IndexedOperator,b::IndexedOperator;h=nothing,kwargs...) = isequal(a.op,b.op)
isscaleequal(a,b;kwargs...) = isequal(a,b)
function has_same(vec1::Vector,vec2::Vector)
    length(vec1) != length(vec2) && return false
    !(issetequal(vec1,vec2)) && return false
    return isequal(counter(vec1),counter(vec2))
end
#function that returns a vector of vector of ops depending on their indices
# so all ops with the same index are returned in the same subvector
function get_ops_of_inds(qmul::QMul; h=nothing, kwargs...)
    inds = get_indices(qmul)
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        filter!(x -> x.specHilb in h,inds)
    end
    dic = Dict{Index,Any}(ind => Any[] for ind in inds)
    allVec = []
    noIndsVec = []
    for op in qmul.args_nc
        if (op isa IndexedOperator && op.ind ∉ inds) || (!(op isa NumberedOperator) && !(op isa IndexedOperator))
            push!(noIndsVec,op)
            continue
        elseif op isa IndexedOperator && op.ind in inds
            push!(dic[op.ind],op)
        end
    end
    for ind in keys(dic)
        push!(allVec,dic[ind])
    end
    push!(allVec,noIndsVec)
    return allVec
end


#function that returns a vector of vector of ops depending on their numbers
# so all ops with the same number are returned in the same subvector
#
# this is still not really correct i think
function get_ops_of_numbs(qmul::QMul; h=nothing, kwargs...)
    numbs = get_numbers(qmul)
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        #filter!(x -> x.specHilb in h,inds)
    end
    dic = Dict{Int64,Any}(numb => Any[] for numb in numbs)
    allVec = []
    noNumbsVec = []
    for op in qmul.args_nc
        if (op isa NumberedOperator && op.numb ∉ numbs) || (!(op isa NumberedOperator) && !(op isa IndexedOperator))
            push!(noNumbsVec,op)
            continue
        elseif op isa NumberedOperator && op.numb in numbs
            push!(dic[op.numb],op)
        end
    end
    for numb in keys(dic)
        push!(allVec,dic[numb])
    end
    push!(allVec,noNumbsVec)
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
            if isscaleequal(avrg,state;kwargs...)
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
