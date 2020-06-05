struct IndexSet{S,T}
    name::S
    lower::T
    upper::T
    function IndexSet{S,T}(name::S, lower::T, upper::T) where {S,T<:Int}
        @assert lower <= upper
        new(name,lower,upper)
    end
end
IndexSet(name::S, lower::T, upper::T) where {S,T} = IndexSet{S,T}(name,lower,upper)
IndexSet(name,upper::Int) = IndexSet(name,1,upper)
IndexSet(name,r::UnitRange{<:Int}) = IndexSet(name,r.start,r.stop)

# Iterator interface
for f in [:IteratorSize, :IteratorEltype, :eltype, :length, :size]
    @eval Base.$(f)(iter::IndexSet, args...) = Base.$(f)(iter.lower:iter.upper, args...)
end

struct Index{I,S} <: Number
    i::I
    set::S
    function Index{I,S}(i::I,set::S) where {I<:Int,S<:IndexSet}
        @assert set.lower <= i <= set.upper
        new(i,set)
    end
end
Index(i::I,set::S) where {I,S} = Index{I,S}(i,set)
Base.isless(i::Index, j::Index) = isless(i.i, j.i)

Base.getindex(I::IndexSet, i::Int) = Index(i,I)
Base.getindex(I::IndexSet, is::Int...) = [Index(i,I) for i in is]

# Iterated values as Index
function Base.iterate(iter::IndexSet, state=(Index(iter.lower, iter), 0))
    element, count = state
    if count >= length(iter)
        return nothing
    end
    if element.i==length(iter)
        return (element, (element, count+1))
    end
    return (element, (Index(element.i + 1, iter), count + 1))
end

default_idx_set() = IndexSet(:DEFAULT, 1, 1)
default_index() = Index(1,default_idx_set())

# _to_symbolic(i::Index) = SymbolicUtils.term(Index, i.i, i.set; type=Int)
# _to_qumulants(I::IndexSet) = I
