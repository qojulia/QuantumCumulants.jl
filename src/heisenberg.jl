"""
    heisenberg(ops::Vector,H::AbstractOperator)
    heisenberg(op::AbstractOperator,H::AbstractOperator)

Compute a set of Heisenberg equations of the operators in `ops`
under the Hamiltonian `H`.
"""
function heisenberg(a::Vector,H; multithread=false)
    if multithread
        lhs = Vector{AbstractOperator}(undef, length(a))
        rhs = Vector{AbstractOperator}(undef, length(a))
        Threads.@threads for i=1:length(a)
            lhs[i] = simplify_operators(a[i])
            rhs[i] = simplify_operators(1.0im*commutator(H,lhs[i];simplify=false))
        end
    else
        lhs = simplify_operators.(a)
        rhs = simplify_operators.([1.0im*commutator(H,a1;simplify=false) for a1=lhs])
    end
    return DifferentialEquation(lhs,rhs)
end
heisenberg(a::AbstractOperator,args...;kwargs...) = heisenberg([a],args...;kwargs...)

"""
    heisenberg(ops::Vector,H::AbstractOperator,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))
    heisenberg(op::AbstractOperator,H::AbstractOperator,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))

Compute the set of equations for the operators in `ops` under the Hamiltonian
`H` and with loss operators contained in `J`. The resulting equation is
equivalent to the Quantum-Langevin equation where noise is neglected.

# Arguments
*`ops::Vector{<:AbstractVector}`: The operators of which the equations are to be computed.
*`H::AbstractOperatr`: The Hamiltonian describing the reversible dynamics of the
    system.
*`J::Vector{<:AbstractOperator}`: A vector containing the collapse operators of
    the system. A term of the form
    ``\\sum_i J_i^\\dagger O J_i - \\frac{1}{2}\\left(J_i^\\dagger J_i O + OJ_i^\\dagger J_i\\right)``
    is added to the Heisenberg equation.

# Optional argumentes
*`Jdagger::Vector=adjoint.(J)`: Vector containing the hermitian conjugates of
    the collapse operators.
*`rates=ones(length(J))`: Decay rates corresponding to the collapse operators in `J`.
"""
function heisenberg(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(length(J)),multithread=false)
    lhs = Vector{AbstractOperator}(undef, length(a))
    rhs = Vector{AbstractOperator}(undef, length(a))
    if multithread
        Threads.@threads for i=1:length(a)
            lhs[i] = simplify_operators(a[i])
            rhs[i] = simplify_operators(1.0im*commutator(H,lhs[i];simplify=false) + _master_lindblad(lhs[i],J,Jdagger,rates))
        end
    else
        for i=1:length(a)
            lhs[i] = simplify_operators(a[i])
            rhs[i] = simplify_operators(1.0im*commutator(H,lhs[i];simplify=false) + _master_lindblad(lhs[i],J,Jdagger,rates))
        end
    end
    return DifferentialEquation(lhs,rhs)
end
function _master_lindblad(a_,J,Jdagger,rates)
    if isa(rates,Vector)
        da_diss = sum(0.5*rates[i]*(Jdagger[i]*commutator(a_,J[i];simplify=false) + commutator(Jdagger[i],a_;simplify=false)*J[i]) for i=1:length(J))
    elseif isa(rates,Matrix)
        da_diss = sum(0.5*rates[i,j]*(Jdagger[i]*commutator(a_,J[j];simplify=false) + commutator(Jdagger[i],a_;simplify=false)*J[j]) for i=1:length(J), j=1:length(J))
    else
        error("Unknown rates type!")
    end
    return simplify_operators(da_diss)
end

function commutator(a::AbstractOperator,b::AbstractOperator; simplify=true, kwargs...)
    # Check on which subspaces each of the operators act
    a_on = acts_on(a)
    b_on = acts_on(b)
    inds = intersect(a_on,b_on)
    isempty(inds) && return zero(a)
    if simplify
        return simplify_operators(a*b - b*a; kwargs...)
    else
        return a*b - b*a
    end
end
