"""
    heisenberg(op::AbstractOperator,H::AbstractOperator)

Compute the Heisenberg equation of the operator `op` under the Hamiltonian `H`.
"""
function heisenberg(a::AbstractOperator,H_)
    H = simplify_operators(H_)
    a_ = simplify_operators(a)
    da = simplify_operators(1.0im*commutator(H,a))
    return DifferentialEquation(a_,da)
end
"""
    heisenberg(ops::Vector,H::AbstractOperator)

Compute a set of Heisenberg equations of the operators in `ops`
under the Hamiltonian `H`.
"""
function heisenberg(a::Vector,H_)
    H = simplify_operators(H_)
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
        da_diss_list = []
        for it=1:length(rates)
            push!(da_diss_list, q_langevin_no_noise(a_,J[it],Jdagger[it],rates[it]))
        end
        da_diss = sum(da_diss_list)
    elseif isa(rates,Matrix)
        error("Matrix not implemented!")
        da_diss = sum(0.5*rates[i,j]*(Jdagger[i]*commutator(a_,J[j]) + commutator(Jdagger[i],a_)*J[j]) for i=1:length(J), j=1:length(J))
    else
        error("Unknown rates type!")
    end
    return apply_comms(da_diss)
end

function q_langevin_no_noise(a_,J_::AbstractOperator,Jdagger_::AbstractOperator,rate)
    c_inds = find_index(rate)[1]
    op_inds = find_index(J_)[2]
    if isempty(c_inds) #is there no IndexedOperator
       return 0.5*rate*(Jdagger_*commutator(a_,J_) + commutator(Jdagger_,a_)*J_)
   elseif length(c_inds) == 1
       if isempty(op_inds)
           error("Need IndexedOperator")
       end
       return Sum(0.5*rate*(Jdagger_*commutator(a_,J_) + commutator(Jdagger_,a_)*J_), op_inds[1])
   elseif length(c_inds) == 2
       if isempty(op_inds)
           error("Need IndexedOperator")
       end
       idx1 = c_inds[1]
       idx2 = c_inds[2]
       J_idx2 = swap_index(J_, op_inds[1], idx2)
       Jd_idx1 = swap_index(Jdagger_, op_inds[1], idx1)
       return Sum(Sum( 0.5*rate*(Jd_idx1*commutator(a_,J_idx2) + commutator(Jd_idx1,a_)*J_idx2),idx1),idx2)
   elseif length(c_inds) >= 3
       error("triple sum not implemented")
   end
end

commutator(a::AbstractOperator,b::AbstractOperator) = apply_comms(a*b - b*a)
function commutator(a::TensorProd,b::TensorProd)
    # Check on which subspaces each of the operators act
    a_on = acts_on(simplify_operators(a))
    b_on = acts_on(simplify_operators(b))
    inds = intersect(a_on,b_on)
    isempty(inds) && return zero(a) # Do not act on same hilbert space
    return apply_comms(a*b - b*a)

    # # Reduce operators to relevant space
    # TODO: optimize the following so one actually gains something
    # red_ab = ⊗([a.args[i]*b.args[i] for i=inds]...)
    # red_ba = ⊗([b.args[i]*a.args[i] for i=inds]...)
    # red_op = apply_comms(red_ab - red_ba)
    # return _expand_operator(red_op,a,b,a_on,b_on,inds,n)
end
# commutator(a::Add,b::TensorProd) = simplify_operators(sum(commutator(a1,b) for a1=a.args))
# commutator(a::TensorProd,b::Add) = simplify_operators(sum(commutator(a,b1) for b1=b.args))
# commutator(a::Add,b::Add) = simplify_operators(sum(commutator(a1,b1) for a1=a.args, b1=b.args))

"""
    acts_on(::AbstractOperator)

Computes on which subspaces the given operator acts nontrivially. Returns a
vector containing the corresponding indices.
"""
acts_on(::BasicOperator) = Int[1]
acts_on(::Identity) = Int[]
acts_on(::Zero) = Int[]
acts_on(::Number) = Int[]
function acts_on(a::Prod)
    args_ = isa(a.args[1],Number) ? a.args[2:end] : a.args
    for a1=args_
        isone(a1) || return Int[1]
    end
    return Int[]
end
function acts_on(a::TensorProd)
    act = Int[]
    for i=1:length(a.args)
        aon = acts_on(a.args[i])
        isempty(aon) || push!(act,i)
    end
    return act
end
function acts_on(a::Add)
    act = collect(Iterators.flatten(acts_on.(a.args)))
    unique!(act)
    return act
end

replace_commutator(a::AbstractOperator,b::AbstractOperator) = (false,a)

apply_comms(a::AbstractOperator) = a
function apply_comms(a::Prod)
    a_ = simplify_operators(a)
    (isone(a_) || iszero(a_)) && return a_

    args = []
    i = 1
    while i < length(a_.args)
        # If we encounter a fundamental operator, check if it matches a canonical commutator
        if isa(a_.args[i], BasicOperator)
            check, arg_ = replace_commutator(a_.args[i],a_.args[i+1])
            push!(args,arg_)
            i += 1 + Int(check)
        elseif isa(a_.args[i], Expression)
            push!(args, apply_comms(a_.args[i]))
            i += 1
        else
            push!(args,a_.args[i])
            i += 1
        end
    end

    # If there was no replacement at the last step, need to add last arg
    if i == length(a_.args)
        push!(args, a_.args[end])
    end

    # TODO: cleaner solution for the recursion here
    op_ = simplify_operators(prod(args))
    if op_ == a # No simplification occurred
        return op_
    else
        return apply_comms(op_)
    end
end
function apply_comms(a::Union{TensorProd,Add})
    a_ = simplify_operators(a)
    if isa(a_, Expression)
        if isa(a_, SumType)
            return simplify_operators(Sum([apply_comms(a1) for a1=a_.args]...,a_.f.index))
        else
            return simplify_operators(a_.f([apply_comms(a1) for a1=a_.args]...))
        end
    else
        return a_
    end
end
apply_comms(a::Number) = a

function simplify_operators(a::Prod)
    # args = eltype(a.args)[]
    args = []
    fac = 1
    for i=1:length(a.args)
        arg_ = simplify_operators(a.args[i])
        iszero(arg_) && return zero(a)
        if isa(arg_, AbstractOperator)
            iszero(arg_) && return arg_
            isone(arg_) || push!(args, arg_)
        else
            fac *= arg_
        end
    end
    fac = isa(fac, SymPy.Sym) ? SymPy.expand(fac) : fac
    iszero(fac) && return Zero()
    isempty(args) && return fac*one(a)

    # Combine pairs where possible
    i = 1
    args2 = AbstractOperator[]
    while i < length(args)
        can_combine, arg_ = combine_prod(args[i],args[i+1])
        # _arg_ = simplify_operators(arg_)
        iszero(arg_) && return zero(a)
        isone(arg_) || push!(args2, arg_)
        i += 1 + Int(can_combine)
    end

    # If there was no replacement at the last step, need to add last arg
    if i == length(args)
        push!(args2, args[end])
    end

    p = prod(args2)
    # TODO: do we need recursion here?
    can_simplify = (a!=fac*p)
    if can_simplify
        return simplify_operators(fac*p)
    else
    isone(fac) && return p
    return fac*p
    end
end
function simplify_operators(a::TensorProd)
    args = []
    fac = 1
    for i=1:length(a.args)
        arg_ = simplify_operators(a.args[i])
        iszero(arg_) && return zero(a)
        if isa(arg_, Prod) && isa(arg_.args[1],Number)
            fac *= arg_.args[1]
            push!(args, prod(arg_.args[2:end]))
        else
            push!(args, arg_)
        end
    end
    out = ⊗(args...)
    # TODO: do we need recursion here?
    fac = isa(fac, SymPy.Sym) ? SymPy.expand(fac) : fac
    iszero(fac) && return zero(a)
    if fac*out == a
        isone(fac) && return out
    return fac*out
    else
        return simplify_operators(fac*out)
    end
end
function simplify_operators(a::Add)
    args = []

    # Simplify each term and skip over zeros
    for i=1:length(a.args)
        arg_ = simplify_operators(a.args[i])
        iszero(arg_) || push!(args, arg_)
    end
    isempty(args) && return zero(a)

    # Collect equal terms
    args2 = []
    skip_index = Int[]
    for i=1:length(args)
        i ∈ skip_index && continue
        arg_ = copy(args[i])
        for j=i+1:length(args)
            can_combine, arg_ = combine_add(arg_,args[j])
            can_combine && push!(skip_index,j)
        end
        iszero(arg_) || push!(args2,arg_)
    end

    isempty(args2) && return zero(a)
    out = sum(args2)
    # TODO: can we avoid this recursion?
    can_simplify = out!=a
    if can_simplify
        return simplify_operators(out)
    else
        return out
    end
end
simplify_operators(a::Number) = a
simplify_operators(a::AbstractOperator) = a
function simplify_operators(a::Transition)
    if a.i==a.j && a.i==a.GS
        inds_ = filter(!isequal(a.GS), [a.inds...])
        return one(a) - sum(Transition(a.label,i,i,a.inds;GS=a.GS) for i=inds_)
    else
        return a
    end
end


"""
    combine_add(a,b)

If possible, combine `a` and `b` under addition.
"""
function combine_add(a::BasicOperator,b::BasicOperator)
    if a==b
        return true, 2*a
    else
        return false, a
    end
end
function combine_add(a::BasicOperator,b::Prod)
    # any([isa(b1,IndexedOperator) for b1=b.args]) && return _combine_add_indexed(a,b)
    if isa(b.args[1],Number) && length(b.args)==2 && a==b.args[2]
        arg_ = copy(b)
        arg_.args[1] += 1
        return true, arg_
    else
        return false, a
    end
end
function combine_add(a::Prod,b::BasicOperator)
    check, arg_ = combine_add(b,a)
    if check
        return true, arg_
    else
        return false, a
    end
end
function combine_add(a::Prod,b::Prod)
    # (any([isa(a1,IndexedOperator) for a1=a.args]) || any([isa(b1,IndexedOperator) for b1=b.args])) && return _combine_add_indexed(a,b)
    if a == b
        return true, 2*a
    elseif isa(a.args[1],Number) && isa(b.args[1],Number) && isequal_prod_args(a.args[2:end],b.args[2:end])
        arg_ = copy(a)
        arg_.args[1] += b.args[1]
        return true, arg_
    elseif isa(a.args[1],Number) && isequal_prod_args(a.args[2:end], b.args)
        arg_ = copy(a)
        arg_.args[1] += 1
        return true, arg_
    elseif isa(b.args[1],Number) && isequal_prod_args(a.args, b.args[2:end])
        arg_ = copy(b)
        arg_.args[1] += 1
        return true, arg_
    else
        return false, a
    end
end
function combine_add(a::TensorProd,b::TensorProd)
    # TODO: Cleaner solution to combine indexed objects
    # _, a_inds = find_index(a)
    # _, b_inds = find_index(b)
    # !(isempty(a_inds) && isempty(b_inds)) && return _combine_add_indexed(a,b)
    check, arg = combine_add(a.args[1],b.args[1])
    if check && (a.args[2:end] == b.args[2:end])
        args = [arg; a.args[2:end]]
        return true, Expression(⊗,args)
    else
        return false, a
    end
end
combine_add(a,b) = (false, a)

combine_prod(a,b) = (false,a)
