"""
Abstract type for equations.
"""
abstract type AbstractEquation{LHS,RHS} end
Base.:(==)(eq1::T,eq2::T) where T<:AbstractEquation = (eq1.lhs==eq2.lhs && eq1.rhs==eq2.rhs)

"""
    DifferentialEquation{LHS,RHS} <: AbstractEquation{LHS,RHS}
    DifferentialEquation(lhs,rhs)

Type defining a differential equation, where `lhs` contains the derivative
and `rhs` the expression.

# Arguments
*`lhs`: The operator or average of which the derivative is taken.
*`rhs`: The expression to which the derivative is equal.

# Fields
*`lhs`: The left-hand side.
*`rhs`: The right-hand side.
"""
mutable struct DifferentialEquation{LHS,RHS} <: AbstractEquation{LHS,RHS}
    lhs::LHS
    rhs::RHS
end

"""
    DifferentialEquationSet{LHS,RHS} <: AbstractEquation{LHS,RHS}
    DifferentialEquationSet(lhs,rhs)

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions.

# Arguments
*`lhs`: Vector of operators or averages of which the derivatives are taken.
*`rhs`: Vector of expressions to which the derivatives are equal.

# Fields
*`lhs`: The left-hand side vector.
*`rhs`: The right-hand side vector.
"""
mutable struct DifferentialEquationSet{LHS,RHS} <: AbstractEquation{LHS,RHS}
    lhs::Vector{LHS}
    rhs::Vector{RHS}
end

Base.getindex(de::DifferentialEquationSet, i::Int) = DifferentialEquation(de.lhs[i],de.rhs[i])
Base.lastindex(de::DifferentialEquationSet) = lastindex(de.lhs)

# Substitution
function substitute(de::DifferentialEquation,dict)
    lhs = substitute(de.lhs, dict)
    rhs = substitute(de.rhs, dict)
    return DifferentialEquation(lhs,rhs)
end
function substitute(de::DifferentialEquationSet,dict)
    lhs = [substitute(l, dict) for l in de.lhs]
    rhs = [substitute(r, dict) for r in de.rhs]
    return DifferentialEquationSet(lhs,rhs)
end

# Simplification
for f in [:simplify_constants,:simplify_operators]
    @eval function $(f)(de::DifferentialEquation;kwargs...)
        lhs = $(f)(de.lhs;kwargs...)
        rhs = $(f)(de.rhs;kwargs...)
        return DifferentialEquation(lhs,rhs)
    end

    @eval function $(f)(de::DifferentialEquationSet;kwargs...)
        lhs = [$(f)(l;kwargs...) for l in de.lhs]
        rhs = [$(f)(r;kwargs...) for r in de.rhs]
        return DifferentialEquationSet(lhs,rhs)
    end
end
