struct Index{S,L,U} <: SymbolicNumber
    name::S
    lower::L
    upper::U
end
Base.hash(i::Index, h::UInt) = hash(i.upper, hash(i.lower, hash(i.name, h)))
Base.isless(i::Index, j::Index) = isless(hash(j), hash(i))
# Base.in(i::Index,j::Index) = isequal(i, j)

default_index() = Index(:DEFAULT, 1, 1)
