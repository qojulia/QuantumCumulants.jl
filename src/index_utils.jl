# get_indices functions
get_indices(x::AvgSums) = get_indices(arguments(x))
function get_indices(term::SymbolicUtils.Term{AvgSym, Nothing})
    get_indices(arguments(term)[1])
end
function get_indices(term::QMul)
    args_nc = copy(term.args_nc)
    filter!(x -> x isa IndexedOperator, args_nc)
    indices = unique(get_indices(args_nc))
end
get_indices(a::IndexedOperator) = [a.ind]
get_indices(vec::Vector) = vcat(get_indices.(vec)...)
get_indices(a::SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}) = unique([a.metadata.ind1,a.metadata.ind2])
get_indices(a::SymbolicUtils.Sym{Parameter,IndexedVariable}) = [a.metadata.ind]
get_indices(x::Number) = []
get_indices(term) = istree(term) ? get_indices(arguments(term)) : []
const Sums = Union{IndexedSingleSum,IndexedDoubleSum}
get_indices(x::Sums) = get_indices(arguments(x))

#Usability functions:
Σ(a,b) = IndexedDoubleSum(a,b)  #Double-Sum here, because if variable a is not a single sum it will create a single sum anyway
Σ(a,b,c;kwargs...) = IndexedDoubleSum(a,b,c;kwargs...)
∑(args...; kwargs...) = Σ(args...; kwargs...)

IndexedOperator(x::indexable,numb::Int64) = NumberedOperator(x,numb)
IndexedVariable(x,numb::Int64) = SingleNumberedVariable(x,numb)
IndexedVariable(x,num1::Int64,num2::Int64;kwargs...) = DoubleNumberedVariable(x,num1,num2;kwargs...)
IndexedVariable(name::Symbol,ind1::Index,ind2::Index;kwargs...) = DoubleIndexedVariable(name,ind1,ind2;kwargs...)

#Numeric Conversion of NumberedOperators
#this only works, if only one of the hilbertspaces is indexed
#function to_numeric(op::NumberedOperator,b::QuantumOpticsBase.CompositeBasis; kwargs...)
#    #keep in mind: this function does not have a check for hilbert-spaces, meaning it can produce wrong output
#    # or error when getting the right basis of b
#    aon = getNumber(op)[1] - 1
#    op_ = _to_numeric(op.op,b.bases[aon];kwargs...)
#    return QuantumOpticsBase.embed(b,aon,op_)
#end
#one needs to call this function, when there are multiple indexed hilbertspaces
#TODO: review and simplify
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
