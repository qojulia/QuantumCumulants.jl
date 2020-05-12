# Conversion to SymbolicUtils
_to_symbolic(op::BasicOperator) = OPERATORS_TO_SYMS[op]
_to_symbolic(t::T) where T<:OperatorTerm = SymbolicUtils.Term{AbstractOperator}(t.f, _to_symbolic.(t.arguments))
_to_symbolic(x::Number) = x
_to_symbolic(x::SymbolicUtils.Symbolic) = x

_to_operator(s::SymbolicUtils.Sym{T}) where T<:AbstractOperator = SYMS_TO_OPERATORS[s]
_to_operator(t::SymbolicUtils.Term) = OperatorTerm(t.f, _to_operator.(t.arguments))
_to_operator(x::Number) = x
_to_operator(s::SymbolicUtils.Symbolic{<:Number}) = s

# Interfacing with SymbolicUtils
SymbolicUtils.promote_symtype(f, Ts::Type{<:AbstractOperator}...) = promote_type(Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:AbstractOperator}, Ts...) = T
# SymbolicUtils._isone(x::SymbolicUtils.Sym{T}) where T<:AbstractOperator = T <: Identity
SymbolicUtils._iszero(x::SymbolicUtils.Sym{T}) where T<:AbstractOperator = T <: Zero

# Symbolic type promotion
for f in [+,-,*,/,^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:Number}) = T
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Number},
                   S::Type{<:AbstractOperator}) = S
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:AbstractOperator}) = promote_type(T,S)
end

Base.one(x::SymbolicUtils.Sym{T}) where T<:BasicOperator = _to_symbolic(one(_to_operator(x)))

function simplify_operators(op::AbstractOperator; rules=SIMPLIFY_OPERATOR_RULES, kwargs...)
    s = _to_symbolic(op)
    s_ = SymbolicUtils.simplify(s; rules=rules, kwargs...)
    (SymbolicUtils.symtype(s_) == Any) && @warn "SymbolicUtils.simplify returned symtype Any; recursion failed!"
    # s_c = SymbolicUtils.simplify(s_, rules=SIMPLIFY_COMMUTATOR_RULES; kwargs...)
    # return s_
    return _to_operator(s_)
end
