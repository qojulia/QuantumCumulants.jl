"""
Abstract type for equations.
"""
abstract type AbstractEquation{LHS,RHS} end
Base.isequal(eq1::T,eq2::T) where T<:AbstractEquation = isequal(eq1.lhs, eq2.lhs) && isequal(eq1.rhs, eq2.rhs)

"""
    DifferentialEquation{LHS,RHS} <: AbstractEquation{LHS,RHS}
    DifferentialEquation(lhs,rhs)

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions.

# Arguments
*`lhs`: Vector of operators or averages of which the derivatives are taken.
*`rhs`: Vector of expressions to which the derivatives are equal.

# Fields
*`lhs`: The left-hand side vector.
*`rhs`: The right-hand side vector.
"""
mutable struct DifferentialEquation{LHS,RHS} <: AbstractEquation{LHS,RHS}
    lhs::Vector{LHS}
    rhs::Vector{RHS}
end

Base.getindex(de::DifferentialEquation, i::Int) = DifferentialEquation([de.lhs[i]],[de.rhs[i]])
Base.lastindex(de::DifferentialEquation) = lastindex(de.lhs)
Base.length(de::DifferentialEquation) = length(de.lhs)

# Substitution
function substitute(de::DifferentialEquation,dict)
    lhs = [substitute(l, dict) for l in de.lhs]
    rhs = [substitute(r, dict) for r in de.rhs]
    return DifferentialEquation(lhs,rhs)
end

# Simplification
for f in [:simplify_constants,:simplify_operators]
    @eval function $(f)(de::DifferentialEquation;kwargs...)
        lhs = [$(f)(l;kwargs...) for l in de.lhs]
        rhs = [$(f)(r;kwargs...) for r in de.rhs]
        return DifferentialEquation(lhs,rhs)
    end
end
