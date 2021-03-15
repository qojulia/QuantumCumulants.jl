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
    @assert length(rates)==length(J)==length(Jdagger)
    J_, Jdagger_, rates_ = _expand_clusters(J,Jdagger,rates)
    if simplify_input
        lhs = map(qsimplify, a)
    else
        lhs = a
    end
    rhs = Vector{QNumber}(undef, length(a))
    if multithread
        Threads.@threads for i=1:length(a)
            rhs_ = 1.0im*commutator(H,lhs[i];simplify=false)
            if !isempty(J)
                rhs_ = rhs_ + _master_lindblad(lhs[i],J_,Jdagger_,rates_)
            end
            rhs[i] = qsimplify(rhs_)
        end
    else
        for i=1:length(a)
            rhs_ = 1.0im*commutator(H,lhs[i];simplify=false)
            if !isempty(J)
                rhs_ = rhs_ + _master_lindblad(lhs[i],J_,Jdagger_,rates_)
            end
            rhs[i] = qsimplify(rhs_)
        end
    end
    if has_cluster(H)
        return scale(HeisenbergEquation(lhs,rhs,H,J_,rates_))
    else
        return HeisenbergEquation(lhs,rhs,H,J_,rates_)
    end
end
heisenberg(a::QNumber,args...;kwargs...) = heisenberg([a],args...;kwargs...)
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
function commutator(a::QNumber,b::QNumber; simplify=true, kwargs...)
    # Check on which subspaces each of the operators act
    a_on = acts_on(a)
    b_on = acts_on(b)
    inds = intersect(a_on,b_on)
    isempty(inds) && return zero(a)
    if simplify
        return qsimplify(a*b + -1*b*a; kwargs...)
    else
        return a*b + -1*b*a
    end
end

# Specialized methods for addition using linearity
function commutator(a::QTerm{<:typeof(+)},b::QNumber; simplify=true, kwargs...)
    args = Any[]
    for arg in a.arguments
        c = commutator(arg,b; simplify=simplify, kwargs...)
        iszero(c) || push!(args, c)
    end
    isempty(args) && return zero(a)
    out = +(args...)
    if simplify
        return qsimplify(out; kwargs...)
    else
        return out
    end
end
function commutator(a::QNumber,b::QTerm{<:typeof(+)}; simplify=true, kwargs...)
    args = Any[]
    for arg in b.arguments
        c = commutator(a,arg; simplify=simplify, kwargs...)
        iszero(c) || push!(args, c)
    end
    isempty(args) && return zero(a)
    out = +(args...)
    if simplify
        return qsimplify(out; kwargs...)
    else
        return out
    end
end
function commutator(a::QTerm{<:typeof(+)},b::QTerm{<:typeof(+)}; simplify=true, kwargs...)
    args = Any[]
    for a_arg in a.arguments
        for b_arg in b.arguments
            c = commutator(a_arg,b_arg; simplify=simplify, kwargs...)
            iszero(c) || push!(args, c)
        end
    end
    isempty(args) && return zero(a)
    out = +(args...)
    if simplify
        return qsimplify(out; kwargs...)
    else
        return out
    end
end

commutator(::Union{T,SymbolicUtils.Symbolic{T}},::QNumber;kwargs...) where T<:Number = 0
commutator(::QNumber,::Union{T,SymbolicUtils.Symbolic{T}};kwargs...) where T<:Number = 0
commutator(::Union{T,SymbolicUtils.Symbolic{T}},::Union{S,SymbolicUtils.Symbolic{S}};kwargs...) where {T<:Number,S<:Number} = 0

function _expand_clusters(J,Jdagger,rates)
    J_ = []
    Jdagger_ = []
    rates_ = []
    for i=1:length(J)
        if J[i] isa Vector
            h = hilbert(J[i][1])
            aon_ = acts_on(J[i][1])
            aon = if aon_ isa Vector
                @assert length(aon_)==1
                aon_[1]
            else
                aon_
            end
            @assert has_cluster(h, aon)
            append!(J_, J[i])
            append!(Jdagger_, Jdagger[i])
            order = h.spaces[get_i(aon)].order
            append!(rates_, [rates[i] for k=1:order])
        else
            push!(J_, J[i])
            push!(Jdagger_, Jdagger[i])
            push!(rates_, rates[i])
        end
    end
    return J_,Jdagger_,rates_
end
