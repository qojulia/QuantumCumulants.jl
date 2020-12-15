"""
Abstract type for equations.
"""
abstract type AbstractEquation{LHS,RHS} end
Base.:(==)(eq1::T,eq2::T) where T<:AbstractEquation = (eq1.lhs==eq2.lhs && eq1.rhs==eq2.rhs)

"""
    DifferentialEquation{LHS,RHS,H,J,R} <: AbstractEquation{LHS,RHS}
    DifferentialEquation(lhs,rhs,H,J,rates)

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions. In addition, it keeps track
of the Hamiltonian, the collapse operators and the corresponding decay rates of
the system.

# Fields
*`lhs`: Vector of operators or averages of which the derivatives are taken.
*`rhs`: Vector of expressions to which the derivatives are equal.
*`hamiltonian`: Operator defining the system Hamiltonian.
*`jumps`: Vector of operators specifying the decay processes.
*`rates`: Decay rates corresponding to the `jumps`.

"""
mutable struct DifferentialEquation{LHS,RHS,H,J,R} <: AbstractEquation{LHS,RHS}
    lhs::Vector{LHS}
    rhs::Vector{RHS}
    hamiltonian::H
    jumps::J
    rates::R
end
Base.hash(eq::DifferentialEquation, h::UInt) = hash(eq.rates, hash(eq.jumps, hash(eq.hamiltonian, hash(eq.rhs, hash(eq.lhs, h)))))
Base.:(==)(eq1::DifferentialEquation,eq2::DifferentialEquation) = hash(eq1)==hash(eq2)

Base.getindex(de::DifferentialEquation, i::Int) = DifferentialEquation([de.lhs[i]],[de.rhs[i]],de.hamiltonian,de.jumps,de.rates)
Base.getindex(de::DifferentialEquation, i_ls::Vector{Int}) = DifferentialEquation([de.lhs[i] for i in i_ls],[de.rhs[i] for i in i_ls],de.hamiltonian,de.jumps,de.rates)
Base.getindex(de::DifferentialEquation, i_ls::UnitRange) = DifferentialEquation([de.lhs[i] for i in i_ls],[de.rhs[i] for i in i_ls],de.hamiltonian,de.jumps,de.rates)
Base.lastindex(de::DifferentialEquation) = lastindex(de.lhs)
Base.length(de::DifferentialEquation) = length(de.lhs)

# Substitution
function substitute(de::DifferentialEquation,dict)
    lhs = [substitute(l, dict) for l in de.lhs]
    rhs = [substitute(r, dict) for r in de.rhs]
    return DifferentialEquation(lhs,rhs,de.hamiltonian,de.jumps,de.rates)
end

# Simplification
for f in [:simplify_constants,:simplify_operators]
    @eval function $(f)(de::DifferentialEquation;kwargs...)
        lhs = [$(f)(l;kwargs...) for l in de.lhs]
        rhs = [$(f)(r;kwargs...) for r in de.rhs]
        return DifferentialEquation(lhs,rhs,de.hamiltonian,de.jumps,de.rates)
    end
end

mutable struct ScaleDifferentialEquation{LHS,RHS,H,J,R,N,ID,INT,D} <: AbstractEquation{LHS,RHS}
    lhs::Vector{LHS}
    rhs::Vector{RHS}
    hamiltonian::H
    jumps::J
    rates::R
    factors::N
    identicals::ID
    interactions::INT
    dictionaries::D
end

Base.hash(eq::ScaleDifferentialEquation, h::UInt) = hash(hash(eq.dictionaries),(hash(eq.interactions),(hash(eq.identicals),(hash(eq.factors),hash(eq.rates, hash(eq.jumps, hash(eq.hamiltonian, hash(eq.rhs, hash(eq.lhs, h)))))))))
Base.:(==)(eq1::ScaleDifferentialEquation,eq2::ScaleDifferentialEquation) = hash(eq1)==hash(eq2)

get_DiffEq_from_scale(de::ScaleDifferentialEquation) = DifferentialEquation(de.lhs,de.rhs,de.hamiltonian,de.jumps,de.rates)

Base.getindex(de::ScaleDifferentialEquation, i::Int) = ScaleDifferentialEquation([de.lhs[i]],[de.rhs[i]],de.hamiltonian,de.jumps,de.rates,de.factors,de.identicals,de.interactions,de.dictionaries)
Base.getindex(de::ScaleDifferentialEquation, i_ls::Vector{Int}) = ScaleDifferentialEquation([de.lhs[i] for i in i_ls],[de.rhs[i] for i in i_ls],de.hamiltonian,de.jumps,de.rates,de.factors,de.identicals,de.interactions,de.dictionaries)
Base.getindex(de::ScaleDifferentialEquation, i_ls::UnitRange) = ScaleDifferentialEquation([de.lhs[i] for i in i_ls],[de.rhs[i] for i in i_ls],de.hamiltonian,de.jumps,de.rates,de.factors,de.identicals,de.interactions,de.dictionaries)
Base.lastindex(de::ScaleDifferentialEquation) = lastindex(de.lhs)
Base.length(de::ScaleDifferentialEquation) = length(de.lhs)

# Substitution
function substitute(de::ScaleDifferentialEquation,dict)
    lhs = [substitute(l, dict) for l in de.lhs]
    rhs = [substitute(r, dict) for r in de.rhs]
    return ScaleDifferentialEquation(lhs,rhs,de.hamiltonian,de.jumps,de.rates,de.factors,de.identicals,de.interactions,de.dictionaries)
end

# Simplification
for f in [:simplify_constants,:simplify_operators]
    @eval function $(f)(de::ScaleDifferentialEquation;kwargs...)
        lhs = [$(f)(l;kwargs...) for l in de.lhs]
        rhs = [$(f)(r;kwargs...) for r in de.rhs]
        return ScaleDifferentialEquation(lhs,rhs,de.hamiltonian,de.jumps,de.rates,de.factors,de.identicals,de.interactions,de.dictionaries)
    end
end
