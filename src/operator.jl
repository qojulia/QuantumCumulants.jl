import Base: +,-,*,==,copy
import LinearAlgebra: ishermitian

"""
    AbstractOperator

Abstract type for all operators and expressions.
"""
abstract type AbstractOperator end# <: Number end

"""
    BasicOperator <: AbstractOperator

Abstract type which is supertype to all fundamental operators, i.e. all operators
that are not `<:Expression`.
"""
abstract type BasicOperator <: AbstractOperator end

Base.iszero(::AbstractOperator) = false
Base.isapprox(a::AbstractOperator,b::AbstractOperator) = isequal(a,b)

"""
    Identity <: BasicOperator
    Identity()
    Base.one(::AbstractOperator)

The identity operator (shown as `𝟙`). This operator leaves any other operator
invariant under multiplication.

# Fields:
*`label`: Symbolic label.
*`id`: Identifier.
"""
mutable struct Identity{L,I} <: BasicOperator
    label::L
    id::I
end
Identity() = Identity(:id,1)
==(::Identity,::Identity) = true

*(a::Identity,b::BasicOperator) = b
*(a::BasicOperator,b::Identity) = a
*(a::Identity,b::Identity) = a
Base.adjoint(a::Identity) = a
Base.one(::BasicOperator) = Identity()
copy(::Identity) = Identity()

"""
    Zero <: BasicOperator
    Zero()
    Base.zero(::AbstractOperator)

The zero operator (shown as `0`). This operator leaves any other operator invariant
under addition. When multiplying any operator with it, however, `Zero()` is returned.
"""
mutable struct Zero{L,I} <: BasicOperator
    label::L
    id::I
end
Zero() = Zero(:Zr,0)
==(::Zero,::Zero) = true

*(a::Zero,b::BasicOperator) = a
*(a::BasicOperator,b::Zero) = b
*(a::Zero,::Identity) = a
*(::Identity,a::Zero) = a
+(a::Zero,b::BasicOperator) = b
+(a::BasicOperator,b::Zero) = a
+(a::Zero,::Zero) = a
Base.adjoint(a::Zero) = a
Base.zero(::BasicOperator) = Zero()
Base.iszero(::Zero) = true
copy(::Zero) = Zero()

*(a::AbstractOperator) = a

ishermitian(a::AbstractOperator) = (a==a')
ishermitian(::Identity) = true
ishermitian(::Zero) = true
