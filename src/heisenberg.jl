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
function heisenberg(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(Int,length(J)),
                    multithread=false,
                    simplify=true,
                    expand=false,
                    order=nothing,
                    mix_choice=maximum,
                    iv=SymbolicUtils.Sym{Real}(:t))

    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    imH = im*H
    if multithread
        Threads.@threads for i=1:length(a)
            rhs_ = commutator(imH,a[i])
            rhs_diss = _master_lindblad(a[i],J,Jdagger,rates)
            rhs[i] = rhs_ + rhs_diss
        end
    else
        for i=1:length(a)
            rhs_ = commutator(imH,a[i])
            rhs_diss = _master_lindblad(a[i],J,Jdagger,rates)
            rhs[i] = rhs_ + rhs_diss
        end
    end

    # Average
    vs = map(average, a)
    rhs_avg = map(average, rhs)
    if simplify
        rhs_avg = map(SymbolicUtils.simplify, rhs_avg)
    end
    rhs = map(undo_average, rhs_avg)

    if expand
        order===nothing && error("Need given order for cumulant expansion!")
        rhs_avg = [cumulant_expansion(r, order; simplify=simplify, mix_choice=mix_choice) for r∈rhs_avg]
    end

    eqs_avg = [Symbolics.Equation(l,r) for (l,r)=zip(vs,rhs_avg)]
    eqs = [Symbolics.Equation(l,r) for (l,r)=zip(a,rhs)]

    varmap = make_varmap(vs, iv)

    return HeisenbergEquation(eqs_avg,eqs,vs,a,H,J,rates,iv,varmap,order)
end
heisenberg(a::QNumber,args...;kwargs...) = heisenberg([a],args...;kwargs...)
heisenberg(a::Vector,H;kwargs...) = heisenberg(a,H,[];Jdagger=[],kwargs...)

function _master_lindblad(a_,J,Jdagger,rates)
    args = Any[]
    if isa(rates,Vector)
        for i=1:length(J)
            c1 = 0.5*rates[i]*Jdagger[i]*commutator(a_,J[i])
            c2 = 0.5*rates[i]*commutator(Jdagger[i],a_)*J[i]
            push_or_append_nz_args!(args, c1)
            push_or_append_nz_args!(args, c2)
        end
    elseif isa(rates,Matrix)
        for i=1:length(J), j=1:length(J)
            c1 = 0.5*rates[i,j]*Jdagger[i]*commutator(a_,J[j])
            c2 = 0.5*rates[i,j]*commutator(Jdagger[i],a_)*J[j]
            push_or_append_nz_args!(args, c1)
            push_or_append_nz_args!(args, c2)
        end
    else
        error("Unknown rates type!")
    end
    isempty(args) && return 0
    return QAdd(args)
end


## Commutator methods
"""
    commutator(a,b)

Computes the commutator `a*b - b*a`.
"""
commutator(a,b) = _commutator(a,b)
_commutator(a, b) = a*b - b*a
commutator(a::QNumber,b::SNuN) = 0
commutator(a::SNuN,b::QNumber) = 0
commutator(::SNuN,::SNuN) = 0
function commutator(a::QSym,b::QSym)
    acts_on(a)==acts_on(b) || return 0
    isequal(a,b) && return 0
    return _commutator(a,b)
end
function commutator(a::QMul,b::QSym)
    aon = acts_on(b)
    idx = findfirst(x->isequal(acts_on(x),aon),a.args_nc)
    idx===nothing && return 0
    return _commutator(a,b)
end
function commutator(a::QSym,b::QMul)
    aon = acts_on(a)
    idx = findfirst(x->isequal(acts_on(x),aon),b.args_nc)
    idx===nothing && return 0
    return _commutator(a,b)
end
function commutator(a::QMul,b::QMul)
    # isequal(a.h, b.h) && return 0
    aon_a = map(acts_on, a.args_nc)
    aon_b = map(acts_on, b.args_nc)
    aon = intersect(aon_a,aon_b)
    isempty(aon) && return 0
    return _commutator(a,b)
end
function commutator(a::QAdd,b::QNumber)
    args = []
    for a_∈a.arguments
        c = commutator(a_,b)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end
function commutator(a::QNumber,b::QAdd)
    args = []
    for b∈b.arguments
        c = commutator(a,b_)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end
function commutator(a::QAdd,b::QAdd)
    args = []
    for a_∈a.arguments, b_∈b.arguments
        c = commutator(a_,b_)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end

function push_or_append_nz_args!(args,c)
    if !SymbolicUtils._iszero(c)
        push!(args, c)
    end
    return args
end
function push_or_append_nz_args!(args,c::QAdd)
    @inbounds for i=1:length(c.arguments)
        push_or_append_nz_args!(args, c.arguments[i])
    end
    return args
end
