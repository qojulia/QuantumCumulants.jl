# getIndices functions
getIndices(x::AvgSums) = getIndices(arguments(x))
function getIndices(term::SymbolicUtils.Term{AvgSym, Nothing})
    arg_ = arguments(term)
    indices = []
    if typeof(arg_[1]) <: QMul
        for arg in arg_[1].args_nc
            if typeof(arg) == IndexedOperator && arg.ind ∉ indices
                push!(indices, arg.ind)
            end
        end
    elseif typeof(arg_[1]) == IndexedOperator
        return [arg_[1].ind]
    end
    return indices
end
getIndices(a::IndexedOperator) = [a.ind]
getIndices(a::SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}) = a.metadata.ind1 == a.metadata.ind2 ? [a.metadata.ind1] : [a.metadata.ind1,a.metadata.ind2]
getIndices(a::SymbolicUtils.Sym{Parameter,IndexedVariable}) = [a.metadata.ind]
getIndices(x::Number) = []
function getIndices(term)
    if istree(term)
        inds = []
        for elem in arguments(term)
            for ind in getIndices(elem)
                if ind ∉ inds
                    push!(inds,ind)
                end
            end
        end
    elseif term isa Vector
        inds = []
        for elem in term
            for ind in getIndices(elem)
                push!(inds,ind)
            end
        end
    else
        return []
    end
    return inds
end
const Sums = Union{IndexedSingleSum,IndexedDoubleSum}
getIndices(x::Sums) = getIndices(arguments(x))


#Usability functions:
Σ(a,b) = IndexedDoubleSum(a,b)  #Double-Sum here, because if variable a is not a single sum it will create a single sum anyway
Σ(a,b,c;kwargs...) = IndexedDoubleSum(a,b,c;kwargs...)
∑(a,b) = Σ(a,b)
∑(a,b,c;kwargs...) = Σ(a,b,c;kwargs...)

IndexedOperator(x::indexable,numb::Int64) = NumberedOperator(x,numb)
function IndexedOperator(x::indexable,num)
    numInt = 0
    try
        numInt = round(Int,num)
    catch
        error("Could not convert given number to Integer!")
    end
    return IndexedOperator(x,numInt)
end
IndexedVariable(x,numb::Int64) = SingleNumberedVariable(x,numb)
function IndexedVariable(x,numb)
    numInt = 0
    try
        numInt = round(Int,numb)
    catch
        error("Could not convert given number to Integer!")
    end
    return IndexedVariable(x,numInt)
end
function IndexedVariable(x,numb1,numb2)
    numInt1 = 0
    numInt2 = 0
    try 
        numInt1 = round(Int,numb1)
        numInt2 = round(Int,numb2)
    catch
        error("Could not convert given numbers to Integer numbers!")
    end
    return IndexedVariable(x,numInt1,numInt2)
end
IndexedVariable(x,num1::Int64,num2::Int64) = DoubleNumberedVariable(x,num1,num2)
IndexedVariable(name::Symbol,ind1::Index,ind2::Index;kwargs...) = DoubleIndexedVariable(name,ind1,ind2;kwargs...)