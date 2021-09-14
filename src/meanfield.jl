"""
    meanfield(ops::Vector,H::QNumber)
    meanfield(op::QNumber,H::QNumber)

    meanfield(ops::Vector,H::QNumber,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))
    meanfield(op::QNumber,H::QNumber,J::Vector;
            Jdagger::Vector=adjoint.(J),rates=ones(length(J)))

Compute the set of equations for the operators in `ops` under the Hamiltonian
`H` and with loss operators contained in `J`. The resulting equation is
equivalent to the Quantum-Langevin equation where noise is neglected.

# Arguments
*`ops::Vector`: The operators of which the equations are to be computed.
*`H::QNumber`: The Hamiltonian describing the reversible dynamics of the
    system.
*`J::Vector{<:QNumber}`: A vector containing the collapse operators of
    the system. A term of the form
    ``\\sum_i J_i^\\dagger O J_i - \\frac{1}{2}\\left(J_i^\\dagger J_i O + OJ_i^\\dagger J_i\\right)``
    is added to the Heisenberg equation.

# Optional argumentes
*`Jdagger::Vector=adjoint.(J)`: Vector containing the hermitian conjugates of
    the collapse operators.
*`rates=ones(length(J))`: Decay rates corresponding to the collapse operators in `J`.
*`multithread=false`: Specify whether the derivation of equations for all operators in `ops`
    should be multithreaded using `Threads.@threads`.
*`simplify=true`: Specify whether the derived equations should be simplified.
*`order=nothing`: Specify to which `order` a [`cumulant_expansion`](@ref) is performed.
    If `nothing`, this step is skipped.
*`mix_choice=maximum`: If the provided `order` is a `Vector`, `mix_choice` determines
    which `order` to prefer on terms that act on multiple Hilbert spaces.
*`iv=SymbolicUtils.Sym{Real}(:t)`: The independent variable (time parameter) of the system.
"""
function meanfield(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(Int,length(J)),
                    multithread=false,
                    simplify=true,
                    order=nothing,
                    mix_choice=maximum,
                    iv=SymbolicUtils.Sym{Real}(:t))

    if rates isa Matrix
        J = [J]; Jdagger = [Jdagger]; rates = [rates]
    end
    J_, Jdagger_, rates_ = _expand_clusters(J,Jdagger,rates)
    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    imH = im*H
    if multithread
        Threads.@threads for i=1:length(a)
            rhs_ = commutator(imH,a[i])
            rhs_diss = _master_lindblad(a[i],J_,Jdagger_,rates_)
            rhs[i] = rhs_ + rhs_diss
        end
    else
        for i=1:length(a)
            rhs_ = commutator(imH,a[i])
            rhs_diss = _master_lindblad(a[i],J_,Jdagger_,rates_)
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

    if order !== nothing
        rhs_avg = [cumulant_expansion(r, order; simplify=simplify, mix_choice=mix_choice) for r∈rhs_avg]
    end

    eqs_avg = [Symbolics.Equation(l,r) for (l,r)=zip(vs,rhs_avg)]
    eqs = [Symbolics.Equation(l,r) for (l,r)=zip(a,rhs)]
    varmap = make_varmap(vs, iv)

    me = MeanfieldEquations(eqs_avg,eqs,vs,a,H,J_,Jdagger_,rates_,iv,varmap,order)
    if has_cluster(H)
        return scale(me;simplify=simplify,order=order,mix_choice=mix_choice)
    else
        return me
    end
end
meanfield(a::QNumber,args...;kwargs...) = meanfield([a],args...;kwargs...)
meanfield(a::Vector,H;kwargs...) = meanfield(a,H,[];Jdagger=[],kwargs...)

function _master_lindblad(a_,J,Jdagger,rates)
    args = Any[]
    for k=1:length(J)
        if isa(rates[k],SymbolicUtils.Sym) || isa(rates[k],Number)
            c1 = 0.5*rates[k]*Jdagger[k]*commutator(a_,J[k])
            c2 = 0.5*rates[k]*commutator(Jdagger[k],a_)*J[k]
            push_or_append_nz_args!(args, c1)
            push_or_append_nz_args!(args, c2)
        elseif isa(rates[k],Matrix)
            for i=1:length(J[k]), j=1:length(J[k])
                c1 = 0.5*rates[k][i,j]*Jdagger[k][i]*commutator(a_,J[k][j])
                c2 = 0.5*rates[k][i,j]*commutator(Jdagger[k][i],a_)*J[k][j]
                push_or_append_nz_args!(args, c1)
                push_or_append_nz_args!(args, c2)
            end
        else
            error("Unknown rates type!")
        end
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
    for b_∈b.arguments
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

function _expand_clusters(J,Jdagger,rates)
    J_ = []
    Jdagger_ = []
    rates_ = []
    for i=1:length(J)
        if (J[i] isa Vector) && !(rates[i] isa Matrix)
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
