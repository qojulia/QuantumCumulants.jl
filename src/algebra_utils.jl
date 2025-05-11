"""
    find_operators(::HilbertSpace, order; names=nothing)

Find all operators that fully define a system up to the given `order`.
"""
function find_operators(h::HilbertSpace, order::Int; names=nothing, kwargs...)
    if names isa Nothing && (unique(typeof.(h.spaces))!=typeof.(h.spaces))
        alph = 'a':'z'
        names_ = Symbol.(alph[1:length(h.spaces)])
    else
        names_ = names
    end
    fund_ops = fundamental_operators(h;names=names_, kwargs...)
    fund_ops = unique([fund_ops;adjoint.(fund_ops)])
    ops = copy(fund_ops)
    for i=2:order
        ops = [ops;fund_ops]
    end

    all_ops = QNumber[]
    for i=1:order
        for c in combinations(ops, i)
            c_ = prod(reverse(c)) # get normal ordering
            iszero(c_) || push!(all_ops, c_)
        end
    end

    filter!(x->!(x isa QAdd), all_ops)
    unique_ops!(all_ops)
    return all_ops
end
find_operators(op::QNumber,args...) = find_operators(hilbert(op),args...)

"""
    fundamental_operators(::HilbertSpace)

Return all fundamental operators for a given Hilbertspace. For example,
a [`FockSpace`](@ref) only has one fundamental operator, `Destroy`.
"""
function fundamental_operators(h::FockSpace,aon::Int=1;names=nothing)
    name = names isa Nothing ? :a : names[aon]
    a = Destroy(h,name)
    return [a]
end
function fundamental_operators(h::NLevelSpace,aon::Int=1;names=nothing)
    sigmas = Transition[]
    lvls = levels(h)
    name = names isa Nothing ? :σ : names[aon]
    for i=1:length(lvls)
        for j=i:length(lvls)
            (i==j) && lvls[i]==ground_state(h) && continue
            s = Transition(h,name,lvls[i],lvls[j])
            push!(sigmas,s)
        end
    end
    return sigmas
end
function fundamental_operators(h::ProductSpace;kwargs...)
    ops = []
    for i=1:length(h.spaces)
        ops_ = fundamental_operators(h.spaces[i],i;kwargs...)
        ops_ = [embed(h,o,i) for o in ops_]
        append!(ops,ops_)
    end
    return ops
end

for T ∈ [:Destroy,:Create,:Transition]
    @eval function embed(h::ProductSpace,op::($T),i)
        fields = [getfield(op, s) for s∈fieldnames($T) if s≠:metadata]
        fields[1] = h
        fields[end] = i
        return $(T)(fields...)
    end
end

"""
    unique_ops(ops)

For a given list of operators, return only unique ones taking into account
their adjoints.
"""
function unique_ops(ops)
    ops_ = deepcopy(ops)
    unique_ops!(ops_)
    return ops_
end

"""
    unique_ops!(ops)

In-place version of [`unique_ops`](@ref).
"""
function unique_ops!(ops)
    hashes = map(hash, ops)
    hashes′ = map(hash, map(_adjoint, ops))
    seen_hashes = UInt[]
    i = 1
    while i <= length(ops)
        if hashes[i] ∈ seen_hashes || hashes′[i] ∈ seen_hashes
            deleteat!(ops, i)
            deleteat!(hashes, i)
            deleteat!(hashes′, i)
        else
            push!(seen_hashes, hashes[i])
            hashes[i]==hashes′[i] || push!(seen_hashes, hashes′[i])
            i += 1
        end
    end
    return ops
end


# Conversion to numerics

"""
    to_numeric(q::QNumber, b::QuantumOpticsBase.Basis; level_map = nothing)
    to_numeric(q::QNumber, state; level_map = nothing)

Convert a symbolic operator `q` to its equivalent numeric (matrix) form on the
basis `b`. The optional argument `level_map` can be set to a dictionary that
specifies how to map levels of a [`Transition`](@ref) to the ones given
in an `NLevelBasis`. **Note:** If the levels of a transition are symbolic,
setting `level_map` is required.

See also: [`numeric_average`](@ref), [`initial_values`](@ref)

Examples
========

julia> to_numeric(Destroy(FockSpace(:fock), :a), FockBasis(10))
Operator(dim=11x11)
  basis: Fock(cutoff=10)[...]

"""
function to_numeric(op::QSym, b::QuantumOpticsBase.Basis; kwargs...)
    check_basis_match(op.hilbert, b; kwargs...)
    return _to_numeric(op, b; kwargs...)
end

_to_numeric(op::Destroy, b::QuantumOpticsBase.FockBasis; kwargs...) = QuantumOpticsBase.destroy(b)
_to_numeric(op::Create, b::QuantumOpticsBase.FockBasis; kwargs...) = QuantumOpticsBase.create(b)
function _to_numeric(op::Pauli, b::QuantumOpticsBase.SpinBasis; kwargs...)
    (b.spinnumber ≠ 1/2) && error("The SpinBasis needs to be Spin-1/2!")
    axis = op.axis
    if axis == 1 # σx
        QuantumOpticsBase.sigmax(b)
    elseif axis == 2 # σy
        QuantumOpticsBase.sigmay(b)
    elseif axis == 3 # σz
        QuantumOpticsBase.sigmaz(b)
    end
end
function _to_numeric(op::Spin, b::QuantumOpticsBase.SpinBasis; kwargs...)
    axis = op.axis
    if axis == 1 # Sx
        QuantumOpticsBase.sigmax(b)*0.5
    elseif axis == 2 # Sy
        QuantumOpticsBase.sigmay(b)*0.5
    elseif axis == 3 # Sz
        QuantumOpticsBase.sigmaz(b)*0.5
    end
end
function _to_numeric(op::Transition, b::QuantumOpticsBase.NLevelBasis; kwargs...)
    i, j = _convert_levels(op; kwargs...)
    return QuantumOpticsBase.transition(b, i, j)
end

function _convert_levels(op; level_map = nothing)
    i, j = op.i, op.j
    if level_map === nothing
        if (!(i isa Number) || !(j isa Number))
            throw(ArgumentError("Mapping from symbolic levels $(i) and $(j) to NLevelBasis requires kwarg level_map to be set"))
        end
        return op.i, op.j  # assume mapping between integers is just equal
    else
        i = level_map[op.i]
        j = level_map[op.j]
        return i, j
    end
end

check_basis_match(h, b; kwargs...) = throw(ArgumentError("Hilbert space $h and basis $b are incompatible!"))
check_basis_match(::FockSpace, ::QuantumOpticsBase.FockBasis; kwargs...) = nothing
check_basis_match(::PauliSpace, ::QuantumOpticsBase.SpinBasis; kwargs...) = nothing
check_basis_match(::SpinSpace, ::QuantumOpticsBase.SpinBasis; kwargs...) = nothing
function check_basis_match(h::NLevelSpace, b::QuantumOpticsBase.NLevelBasis; kwargs...)
    if length(h.levels) != length(b)
        throw(ArgumentError("Hilbert space $h and basis $b have incompatible levels!"))
    end
end

function check_basis_match(h::ProductSpace, b::QuantumOpticsBase.CompositeBasis; ranges=[], kwargs...)
    if length(h.spaces) != length(b.bases) && (isempty(ranges) || sum(ranges) != length(b.bases))
        throw(ArgumentError("Hilbert space $h and basis $b don't have the same number of subspaces!
             If you use indices, specify the `ranges` kwarg."))
    end
    if isempty(ranges)
        inds = [1:1:length(h.spaces);]
    else
        inds = [sum(ranges[1:i]) for i=1:length(ranges)]
    end
    b_r = [b.bases[i] for i in inds]
    for (h_, b_) ∈ zip(h.spaces, b_r)
        check_basis_match(h_, b_; ranges=ranges)
    end
end


# Composite bases
function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis; kwargs...)
    check_basis_match(op.hilbert, b; kwargs...)
    aon = acts_on(op)
    op_num = _to_numeric(op, b.bases[aon]; kwargs...)
    return QuantumOpticsBase.LazyTensor(b, aon, op_num)
end

# Symbolic expressions
function to_numeric(op::QTerm, b::QuantumOpticsBase.Basis; kwargs...)
    f = SymbolicUtils.operation(op)
    return _to_numeric_term(f, op, b; kwargs...)
end

function _to_numeric_term(f::Function, op, b; kwargs...)
    args = SymbolicUtils.arguments(op)
    return f((to_numeric(arg, b; kwargs...) for arg in args)...)
end

function _to_numeric_term(::typeof(*), op::QTerm, b::QuantumOpticsBase.Basis; kwargs...)
    args = SymbolicUtils.arguments(op)
    factor = 1
    args_num = Any[]
    for arg in args
        if arg isa Number
            factor *= arg
        else
            push!(args_num, to_numeric(arg, b; kwargs...))
        end
    end

    if length(args_num) == 0
        return factor * one(b)
    end

    return *(factor, args_num...)
end

function to_numeric(x::Number, b::QuantumOpticsBase.Basis; kwargs...)
    op = _lazy_one(b)*x
    return op
end
_lazy_one(b::QuantumOpticsBase.Basis) = one(b)
function _lazy_one(b::QuantumOpticsBase.CompositeBasis)
    LazyTensor(b, [1:length(b.bases);], Tuple(one(b_) for b_ in b.bases))
end


"""
    numeric_average(avg::Average, state; level_map = nothing)
    numeric_average(q::QNumber, state; level_map = nothing)

From a symbolic average `avg` or operator `q`, compute the corresponding
numerical average value with the given quantum state `state`. This state
can either be of type `QuantumOpticsBase.StateVector` or `QuantumOpticsBase.Operator`.

See also: [`initial_values`](@ref), [`to_numeric`](@ref)
"""
function numeric_average(avg::Average, state; kwargs...)
    op = undo_average(avg)
    return numeric_average(op, state; kwargs...)
end
to_numeric(op::QNumber, state; kwargs...) = to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)

function numeric_average(op::QNumber, state; kwargs...)
    op_num = to_numeric(op, state; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end

function _conj(v::Average)
    arg = v.arguments[1]
    adj_arg = adjoint(arg)
    if has_cluster(arg)
        aons, N = get_cluster_stuff(hilbert(arg))
        names = get_names(arg)
        return substitute_redundants(_average(adj_arg), aons, names)
    else
        return _average(adj_arg)
    end
end
function _conj(v::T) where T <: SymbolicUtils.Symbolic
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(_conj, SymbolicUtils.arguments(v))
        return SymbolicUtils.maketerm(T, f, args, TermInterface.metadata(v))
    else
        return conj(v)
    end
end
_conj(x::Number) = conj(x)

_adjoint(op::QNumber) = adjoint(op)
_adjoint(s::SymbolicUtils.Symbolic{<:Number}) = _conj(s)
_adjoint(x) = adjoint(x)
