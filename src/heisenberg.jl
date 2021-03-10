"""
    heisenberg(ops::Vector,H::QNumber)
    heisenberg(op::QNumber,H::QNumber)

    heisenberg(ops::Vector,H::QNumber,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))
    heisenberg(op::QNumber,H::QNumber,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))

Compute the set of equations for the operators in `ops` under the Hamiltonian
`H` and with loss operators contained in `J`. The resulting equation is
equivalent to the Quantum-Langevin equation where noise is neglected.

# Arguments
*`ops::Vector{<:AbstractVector}`: The operators of which the equations are to be computed.
*`H::AbstractOperatr`: The Hamiltonian describing the reversible dynamics of the
    system.
*`J::Vector{<:QNumber}`: A vector containing the collapse operators of
    the system. A term of the form
    ``\\sum_i J_i^\\dagger O J_i - \\frac{1}{2}\\left(J_i^\\dagger J_i O + OJ_i^\\dagger J_i\\right)``
    is added to the Heisenberg equation.

# Optional argumentes
*`Jdagger::Vector=adjoint.(J)`: Vector containing the hermitian conjugates of
    the collapse operators.
*`rates=ones(length(J))`: Decay rates corresponding to the collapse operators in `J`.
"""
function heisenberg(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(length(J)),
                    multithread=false,simplify_input=false)
    if simplify_input
        lhs = map(qsimplify, a)
    else
        lhs = a
    end
    rhs = Vector{QSymbolic}(undef, length(a))
    if multithread
        Threads.@threads for i=1:length(a)
            rhs_ = 1.0im*commutator(H,lhs[i];simplify=false)
            if !isempty(J)
                rhs_ = rhs_ + _master_lindblad(lhs[i],J,Jdagger,rates)
            end
            rhs[i] = qsimplify(rhs_)
        end
    else
        for i=1:length(a)
            rhs_ = 1.0im*commutator(H,lhs[i];simplify=false)
            if !isempty(J)
                rhs_ = rhs_ + _master_lindblad(lhs[i],J,Jdagger,rates)
            end
            rhs[i] = qsimplify(rhs_)
        end
    end
    return HeisenbergEquation(lhs,rhs,H,J,rates)
end
heisenberg(a::QSymbolic,args...;kwargs...) = heisenberg([a],args...;kwargs...)
heisenberg(a::Vector,H;kwargs...) = heisenberg(a,H,[];Jdagger=[],kwargs...)

function _master_lindblad(a_,J,Jdagger,rates)
    if isa(rates,Vector)
        da_diss = sum(0.5*rates[i]*(Jdagger[i]*commutator(a_,J[i];simplify=false) + commutator(Jdagger[i],a_;simplify=false)*J[i]) for i=1:length(J))
    elseif isa(rates,Matrix)
        da_diss = sum(0.5*rates[i,j]*(Jdagger[i]*commutator(a_,J[j];simplify=false) + commutator(Jdagger[i],a_;simplify=false)*J[j]) for i=1:length(J), j=1:length(J))
    else
        error("Unknown rates type!")
    end
    return qsimplify(da_diss)
end

"""
    commutator(a,b; simplify=true, kwargs...)

Computes the commutator `a*b - b*a` of `a` and `b`. If `simplify` is `true`, the
result is simplified using the [`qsimplify`](@ref) function. Further
keyword arguments are passed to simplification.
"""
function commutator(a::QSymbolic,b::QSymbolic; simplify=true, kwargs...)
    if SymbolicUtils.istree(a) && SymbolicUtils.operation(a)===(+)
        args = Any[]
        for arg in a.arguments
            c = commutator(arg,b; simplify=simplify, kwargs...)
            iszero(c) || push!(args, c)
        end
        isempty(args) && return zero(a)
        out = +(args...)
    elseif SymbolicUtils.istree(b) && SymbolicUtils.operation(b)===(+)
        args = Any[]
        for arg in b.arguments
            c = commutator(a,arg; simplify=simplify, kwargs...)
            iszero(c) || push!(args, c)
        end
        isempty(args) && return zero(a)
        out = +(args...)
    else
        # Check on which subspaces each of the operators act
        a_on = acts_on(a)
        b_on = acts_on(b)
        inds = intersect(a_on,b_on)
        isempty(inds) && return zero(a)
        out = a*b + -1*b*a
    end
    if simplify
        return qsimplify(out; kwargs...)
    else
        return out
    end
end
