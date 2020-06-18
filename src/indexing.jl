struct Index{S,L,U} <: SymbolicNumber
    name::S
    lower::L
    upper::U
end
Base.hash(i::Index, h::UInt) = hash(i.upper, hash(i.lower, hash(i.name, h)))
Base.isless(i::Index, j::Index) = isless(hash(j), hash(i))
# Base.in(i::Index,j::Index) = isequal(i, j)

default_index() = Index(:DEFAULT, 1, 1)

_to_symbolic(i::Index) = SymbolicUtils.term(index, i.name, i.lower, i.upper; type=Int)

neq_inds_prod(args,inds) = OperatorTerm(neq_inds_prod, [args,inds])
# neq_inds_prod(args) = OperatorTerm(neq_inds_prod, [args,])

_to_symbolic(args::Vector) = _to_symbolic.(args)
_to_qumulants(args::Vector) = _to_qumulants.(args)
