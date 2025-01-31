import Base: sum

# QNumber sums
struct Sum{T<:QNumber,I,M} <: QTerm
    term::T
    index::I
    metadata::M
end

# constructor for nested sums
function Sum(term, indices::Index...)
    # sort indices to make sure we always construct nested sums in the same order
    sorted_indices = sort!([indices...], by=i -> i.name)
    return _nested_sum(term, sorted_indices...)
end

# this method is required to avoid ambiguity even though it's the exact same as above
function Sum(term::QNumber, indices::Index...)
    # sort indices to make sure we always construct nested sums in the same order
    sorted_indices = sort!([indices...], by=i -> i.name)
    return _nested_sum(term, sorted_indices...)
end

function _nested_sum(term, index1::Index, index2::Index, remaining_indices...)
    s_inner = Sum(term, index1)
    return _nested_sum(s_inner, index2, remaining_indices...)
end
_nested_sum(s, index::Index) = Sum(s, index)


hilbert(s::Sum) = hilbert(s.term)

# Summation over CNumbers needs to be a SymbolicUtils.Symbolic{<:CNumber}
struct CSumSym <: CNumber end

csum(args...) = Sum(args...)

# SymbolicUtils.maketerm(::Type{T}, ::typeof(csum), args, metadata) where T = csum(args...)
# SymbolicUtils.maketerm(::Type{<:SymbolicUtils.BasicSymbolic}, ::typeof(csum), args, metadata) = csum(args...)

for f in [:*, :+, :-]
    @eval SymbolicUtils.promote_symtype(::typeof($f), ::Type{CSumSym}, ::Type{CSumSym}) = CNumber
end

function Sum(term::SymbolicUtils.Symbolic{<:Number}, index::Index; metadata = nothing)
    # TODO: don't ignore metadata here
    # TODO: printing for CNumber sums
    if !has_index(term, index)
        return index.range * term
    end

    if SymbolicUtils.is_operation(*)(term)
        has_equality_for_index, to_index = find_equality_for_index(term, index)
        if has_equality_for_index
            return change_index(term, index, to_index)
        end
    end

    return SymbolicUtils.Term{CSumSym}(csum, [term, index])
end


const Σ = Sum
const ∑ = Sum

# Basic methods
Base.isequal(s1::Sum, s2::Sum) = isequal(s1.index, s2.index) && isequal(s1.term, s2.term)
Base.adjoint(s::Sum) = Sum(adjoint(s.term), s.index, s.metadata)

# Symbolics interface
TermInterface.iscall(::Sum) = true
SymbolicUtils.operation(::Sum) = Sum
SymbolicUtils.arguments(s::Sum) = [s.term, s.index]
SymbolicUtils.maketerm(::Type{<:Sum}, ::typeof(Sum), args, metadata) = Sum(args...; metadata)

# has_index
"""
    has_index(expr, i::Index)

Check if an expression has a specific index.

"""
has_index(x, i::Index) = false
has_index(v::IndexedVariable, i::Index) = isequal(v.ind, i)
has_index(v::DoubleIndexedVariable, i::Index) = isequal(v.ind1, i) || isequal(v.ind2, i)
has_index(op::IndexedOperator, i::Index) = isequal(op.ind, i)
has_index(i::Index, j::Index) = isequal(i, j)

has_index(s::Sum, i::Index) = isequal(s.index, i) || has_index(s.term, i)

function has_index(t::SymbolicUtils.Symbolic, i::Index)
    if !TermInterface.iscall(t)
        return false
    end

    return has_index(SymbolicUtils.arguments(t), i)
end
# TODO: specialization for QAdd, QMul, etc.
has_index(t::QTerm, i::Index) = has_index(SymbolicUtils.arguments(t), i)

function has_index(args::Vector, i::Index)
    for arg in args
        if has_index(arg, i)
            return true
        end
    end
    return false
end


function get_indices(x)
    # TODO: can we improve the typing of the set here?
    indices = Set{Union{Index, Int}}()
    get_indices!(indices, x)
end

function get_indices!(indices, x)
    !TermInterface.iscall(x) && return indices

    get_indices!(indices, SymbolicUtils.arguments(x))
    return indices
end

function get_indices!(indices, x::Vector)
    for arg in x
        get_indices!(indices, arg)
    end
    return indices
end

get_indices!(indices, i::Index) = push!(indices, i)
get_indices!(indices, op::IndexedOperator) = push!(indices, op.ind)
get_indices!(indices, eqs::AbstractMeanfieldEquations) = get_indices!(indices, eqs.equations)
function get_indices!(indices, eq::Symbolics.Equation)
    get_indices!(indices, eq.lhs)
    get_indices!(indices, eq.rhs)
    return indices
end

# find_equality_for_index --  checks for occurrence of expressions such as i == j
# which allows to e.g. simplify sums
function find_equality_for_index(t, index::Index)
    if !TermInterface.iscall(t)
        return false, nothing
    end

    args = SymbolicUtils.arguments(t)
    if (TermInterface.operation(t) === (==)) && any(isequal(index), args)
        if length(args) != 2
            throw(error("Equality with more than two arguments encountered! Please report this issue!"))
        end

        element_index = findfirst(!isequal(index), args)
        to_index = args[element_index]
        return true, to_index
    end

    for arg in args
        has_equality_for_index, to_index = find_equality_for_index(arg, index)
        if has_equality_for_index
            return true, to_index
        end
    end

    return false, nothing
end



# Construction of sums with QSyms -- order should be QSym < QMul < QAdd < Sum
Sum(a::QSym, index::Index) = index.range * a

function Sum(a::IndexedOperator, index::Index)
    if !has_index(a, index)
        return index.range * a
    end
    return Sum(a, index, nothing)
end

function Sum(t::T, index::I) where {T<:QMul, I<:Index}
    if !has_index(t, index)
        return index.range * t
    end

    has_equality_for_index, to_index = find_equality_for_index(t, index)
    if has_equality_for_index
        return change_index(t, index, to_index)
    end

    if !has_index(t.args_nc, index)
        # NOTE: this check here is only valid since we resolve i == j above
        # otherwise the condition here is not specific enough to treat e.g. Sum((i == j) * σᵢ, j)
        # appropriately;
        # TODO: make sure this actually works; if it gives us trouble, we can skip this "simplification"
        # here since averaging should then lead to similar simplifications for CNumbers only anyway
        if length(t.args_nc) == 1
            return Sum(t.arg_c, index) * t.args_nc[1]
        else
            return Sum(t.arg_c, index) * *(t.args_nc...)
        end
    end

    return Sum(t, index, nothing)
end

function Sum(s::Sum, index::Index)
    if !has_index(s, index)
        return index.range * s
    end
    return Sum(s, index, nothing)
end

function Sum(t::QAdd, index::Index)
    args = SymbolicUtils.arguments(t)

    # check for delta_ij terms
    for (i, arg) in enumerate(args)
        has_equality_for_index, to_index = find_equality_for_index(arg, index)
        if has_equality_for_index
            new_arg = change_index(arg, index, to_index)
            new_args = vcat(args[1:i-1], new_arg, args[i+1:end])
            return Sum(+(new_args...), index)
        end
    end

    # check if any arguments don't have the summation index
    for (i, arg) in enumerate(args)
        if !has_index(arg, index)
            new_arg = index.range * arg
            remaining_args = vcat(args[1:i-1], args[i+1:end])
            isempty(remaining_args) && return new_arg
            return new_arg + Sum(+(remaining_args...), index)
        end
    end

    return Sum(t, index, nothing)
end




# function Sum(t::QAdd, index::Index)
#     args = [Sum(arg, index) for arg in SymbolicUtils.arguments(t)]
#     if length(args) == 1
#         return args[1]
#     end
#     return +(args...)
# end

# CNumbers
function Sum(t::SymbolicUtils.Symbolic, index::Index)
    if !has_index(t, index)
        return index.range * t
    end
    return Sum(t, index, nothing)
end
Sum(t::Number, index::Index) = index.range * t

# Basic algebra
function +(s1::Sum, s2::Sum)
    # TODO: check for range of indices as well!
    if isequal(s1.index, s2.index)
        term = s1.term + s2.term
        index = s1.index
    elseif !has_index(s2.term, s1.index)
        term = s1.term + change_index(s2.term, s2.index, s1.index)
        index = s1.index
    elseif !has_index(s1.term, s2.index)
        term = change_index(s1.term, s1.index, s2.index)
        index = s2.index
    else
        throw(error("David was lazy"))
    end

    Sum(term, index)
end

# function +(s1::QNumber, s2::Sum)
#     return QAdd([s1, s2])
# end
# function +(s1::Sum, s2::QNumber)
#     return QAdd([s1, s2])
# end

function *(s1::Sum, s2::Sum)
    if isequal(s1.index, s2.index)
        new_index = Index(s2.index.hilbert, gensym(s2.index.name), s2.index.range, s2.index.aon)
        return s1 * change_index(s2, s2.index, new_index)
    end

    return Sum(Sum(s1.term * s2.term, s2.index), s1.index)
end

function *(a::QSym, s::Sum)
    return Sum(a * s.term, s.index)
end
function *(s::Sum, a::QSym)
    return Sum(s.term * a, s.index)
end

function *(a::IndexedOperator, s::Sum)
    if isequal(a.ind, s.index)
        new_index = Index(a.ind.hilb, gensym(a.ind.name), a.ind.range, a.ind.aon)
        return change_index(a, a.ind, new_index) * s
    end
    return Sum(a * s.term, s.index)
end
function *(s::Sum, a::IndexedOperator)
    if isequal(a.ind, s.index)
        new_index = Index(a.ind.hilb, gensym(a.ind.name), a.ind.range, a.ind.aon)
        return s * change_index(a, a.ind, new_index)
    end
    return Sum(s.term * a, s.index)
end

function *(s::Sum, t::QTerm)
    if has_index(t, s.index)
        new_index = Index(s.index.hilb, gensym(s.index.name), s.index.range, s.index.aon)
        return s * change_index(t, s.index, new_index)
    end
    return Sum(s.term * t, s.index)
end
function *(t::QTerm, s::Sum)
    if has_index(t, s.index)
        new_index = Index(s.index.hilb, gensym(s.index.name), s.index.range, s.index.aon)
        return change_index(t, s.index, new_index) * s
    end
    return Sum(t * s.term, s.index)
end

function *(v::SNuN, s::Sum)
    if has_index(v, s.index)
        new_index = Index(s.index.hilb, gensym(s.index.name), s.index.range, s.index.aon)
        return change_index(v, s.index, new_index) * s
    end
    return Sum(v * s.term, s.index)
end
*(s::Sum, v::SNuN) = v * s

# averaging
average(s::Sum) = Sum(average(s.term), s.index)
average(s::Sum, order) = Sum(average(s.term, order), s.index)


# TODO: move this somewhere else where it makes more sense
function _push_lindblad_term!(args::Vector, a::QNumber, rate::SymbolicUtils.Symbolic{<:IndexedParameterSym}, J::QNumber, Jdagger::QNumber)
    indices = get_indices(rate)
    @assert length(indices) == 1
    i = indices[1]
    c1 = 0.5*rate*Jdagger*commutator(a,J)
    c2 = 0.5*rate*commutator(Jdagger,a)*J
    push_or_append_nz_args!(args, Sum(c1, i))
    push_or_append_nz_args!(args, Sum(c2, i))
    return nothing
end


function _push_lindblad_term!(args::Vector, a::QNumber, rate::SymbolicUtils.Symbolic{<:IndexedParameterSym}, J::Vector, Jdagger::Vector)
    c1 = 0.5*rate*Jdagger[1]*commutator(a,J[2])
    c2 = 0.5*rate*commutator(Jdagger[2],a)*J[1]

    indices = get_indices(rate)
    push_or_append_nz_args!(args, Sum(c1, indices...))
    push_or_append_nz_args!(args, Sum(c2, indices...))

    return nothing
end