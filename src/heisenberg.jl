"""
    heisenberg(op::AbstractOperator,H::AbstractOperator)

Compute the Heisenberg equation of the operator `op` under the Hamiltonian `H`.
"""
function heisenberg(a::AbstractOperator,H)
    a_ = simplify_operators(a)
    da = simplify_operators(1.0im*commutator(H,a))
    return DifferentialEquation(a_,da)
end
"""
    heisenberg(ops::Vector,H::AbstractOperator)

Compute a set of Heisenberg equations of the operators in `ops`
under the Hamiltonian `H`.
"""
function heisenberg(a::Vector,H)
    lhs = simplify_operators.(a)
    rhs = simplify_operators.([1.0im*commutator(H,a1) for a1=lhs])
    return DifferentialEquationSet(lhs,rhs)
end

"""
    heisenberg(op::AbstractOperator,H::AbstractOperator,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))

Compute the equation for the operator `op` under the Hamiltonian `H` and with
loss operators contained in `J`. The resulting equation is equivalent to the
Quantum-Langevin equation where noise is neglected.

# Arguments
*`op::AbstractOperator`: The operator of which the equation is to be computed.
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
function heisenberg(a::AbstractOperator,H,J::Vector;Jdagger::Vector=adjoint.(J),rates=ones(length(J)))
    he = heisenberg(a,H)
    a_ = he.lhs
    da_diss = _master_lindblad(a_,J,Jdagger,rates)
    da_ = simplify_operators(he.rhs + da_diss)
    return DifferentialEquation(a_,da_)
end

"""
    heisenberg(ops::Vector,H::AbstractOperator,J::Vector;
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
function heisenberg(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(length(J)))
    he = heisenberg(a,H)
    lhs = he.lhs
    da_diss = [_master_lindblad(a_,J,Jdagger,rates) for a_=lhs]
    da_ = [simplify_operators(he.rhs[i] + da_diss[i]) for i=1:length(lhs)]
    return DifferentialEquationSet(lhs,da_)
end
function _master_lindblad(a_,J,Jdagger,rates)
    if isa(rates,Vector)
        da_diss = sum(0.5*rates[i]*(Jdagger[i]*commutator(a_,J[i]) + commutator(Jdagger[i],a_)*J[i]) for i=1:length(J))
    elseif isa(rates,Matrix)
        da_diss = sum(0.5*rates[i,j]*(Jdagger[i]*commutator(a_,J[j]) + commutator(Jdagger[i],a_)*J[j]) for i=1:length(J), j=1:length(J))
    else
        error("Unknown rates type!")
    end
    return simplify_operators(da_diss)
end

commutator(a::AbstractOperator,b::AbstractOperator) = simplify_operators(a*b - b*a)
function commutator(a::OperatorTerm{<:typeof(⊗)},b::OperatorTerm{<:typeof(⊗)})
    # Check on which subspaces each of the operators act
    a_on = acts_on(simplify_operators(a))
    b_on = acts_on(simplify_operators(b))
    inds = intersect(a_on,b_on)
    isempty(inds) && return 0#zero(a) # Do not act on same hilbert space
    return simplify_operators(a*b - b*a)
end

"""
    acts_on(::AbstractOperator)

Computes on which subspaces the given operator acts nontrivially. Returns a
vector containing the corresponding indices.
"""
acts_on(::BasicOperator) = Int[1]
acts_on(::Identity) = Int[]
acts_on(::Zero) = Int[]
acts_on(x) = Int[]
function acts_on(a::OperatorTerm{<:typeof(*)})
    args = filter(!SymbolicUtils.sym_isa(Number), a.arguments)
    for a1=args
        isa(a1,Identity) || return Int[1]
    end
    return Int[]
end
function acts_on(a::OperatorTerm{<:typeof(⊗)})
    act = Int[]
    for i=1:length(a.arguments)
        aon = acts_on(a.arguments[i])
        isempty(aon) || push!(act,i)
    end
    return act
end
function acts_on(a::OperatorTerm)
    act = collect(Iterators.flatten(acts_on.(a.arguments)))
    unique!(act)
    return act
end
