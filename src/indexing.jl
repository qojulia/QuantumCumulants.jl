const IDX_TO_SYMS = Dict{SymbolicNumber,SymbolicUtils.Sym}()
const SYMS_TO_IDX = Dict{SymbolicUtils.Sym,SymbolicNumber}()

### Indices

struct Index{T<:Int,S,U} <: SymbolicNumber{T}
    name::S
    count::U
    function Index{T,S,U}(name::S,count::U) where {T<:Int,S,U}
        idx = new(name,count)
        if !haskey(IDX_TO_SYMS, idx)
            sym = SymbolicUtils.Sym{Index}(gensym(:Index))
            IDX_TO_SYMS[idx] = sym
            SYMS_TO_IDX[sym] = idx
        end
        return idx
    end
end
Index{T}(name::S,count::U) where {T,S,U} = Index{T,S,U}(name,count)
Index(name,count) = Index{Int}(name,count)
_to_symbolic(idx::Index) = IDX_TO_SYMS[idx]
_to_qumulants(t::SymbolicUtils.Sym{T}) where T<:Index = SYMS_TO_IDX[t]

Base.hash(i::Index, h::UInt) = hash(i.count, hash(i.name, h))
Base.isless(i::Index, j::Index) = isless(hash(j), hash(i))
Base.isequal(i::Index, j::Index) = isequal(hash(j), hash(i))


### Indexed operators

for Tname in [:Destroy,:Create]
    Name = Symbol(:Indexed,Tname)
    @eval struct $(Name){H<:HilbertSpace,S,A,IND} <: BasicOperator
        hilbert::H
        name::S
        aon::A
        index::IND
        function $(Name){H,S,A,IND}(hilbert::H,name::S,aon::A,index::IND) where {H,S,A,IND}
            @assert has_hilbert(FockSpace,hilbert,aon)
            op = new(hilbert,name,aon,index)
            if !haskey(OPERATORS_TO_SYMS, op)
                sym = SymbolicUtils.Sym{$(Name)}(gensym(nameof($Name)))
                OPERATORS_TO_SYMS[op] = sym
                SYMS_TO_OPERATORS[sym] = op
            end
            return op
        end
    end
    @eval $(Name)(hilbert::H,name::S,aon::A,index::IND) where {H,S,A,IND} = $(Name){H,S,A,IND}(hilbert,name,aon,index)
    @eval Base.getindex(a::$Tname,i::Union{Index,Int}) = $(Name)(a.hilbert,a.name,a.aon,i)
    @eval Base.:(==)(a::T, b::T) where T<:($Name) = (a.hilbert==b.hilbert && a.name==b.name && a.aon==b.aon && isequal(a.index, b.index))
    @eval Base.hash(a::$Name, h::UInt) = hash(a.hilbert, hash(a.name, hash(a.aon, hash(a.index, h))))
end
Base.adjoint(a::IndexedDestroy) = IndexedCreate(a.hilbert,a.name,a.aon,a.index)
Base.adjoint(a::IndexedCreate) = IndexedDestroy(a.hilbert,a.name,a.aon,a.index)

struct IndexedTransition{H,S,I,A,IND} <: BasicOperator
    hilbert::H
    name::S
    i::I
    j::I
    aon::A
    index::IND
    function IndexedTransition{H,S,I,A,IND}(hilbert::H,name::S,i::I,j::I,aon::A,index::IND) where {H,S,I,A,IND}
        @assert has_hilbert(NLevelSpace,hilbert,aon)
        @assert i∈levels(hilbert,aon) && j∈levels(hilbert,aon)
        op = new(hilbert,name,i,j,aon,index)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{IndexedTransition}(gensym(:IndexedTransition))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
IndexedTransition(hilbert::H,name::S,i::I,j::I,aon::A,index::IND) where {H,S,I,A,IND} = IndexedTransition{H,S,I,A,IND}(hilbert,name,i,j,aon,index)
Base.getindex(s::Transition,k::Union{Index,Int}) = IndexedTransition(s.hilbert,s.name,s.i,s.j,s.aon,k)
Base.adjoint(s::IndexedTransition) = IndexedTransition(s.hilbert,s.name,s.j,s.i,s.aon,s.index)
Base.:(==)(t1::IndexedTransition,t2::IndexedTransition) = (t1.hilbert==t2.hilbert && t1.name==t2.name && t1.i==t2.i && t1.j==t2.j && isequal(t1.index,t2.index))
Base.hash(t::IndexedTransition, h::UInt) = hash(t.hilbert, hash(t.name, hash(t.i, hash(t.j, hash(t.aon, hash(t.index, h))))))
nip(args::AbstractOperator...) = OperatorTerm(nip, [args...])
nip(args::Vector{<:AbstractOperator}) = OperatorTerm(nip, args)
nip(args::Union{SymbolicUtils.Symbolic,Number}...) = SymbolicUtils.Term{AbstractOperator}(nip, [args...])

### Simplification functions
function commute_bosonic_idx(a::SymbolicUtils.Symbolic,b::SymbolicUtils.Symbolic)
    return _to_symbolic(commute_bosonic_idx(_to_qumulants(a), _to_qumulants(b)))
end
function commute_bosonic_idx(a,b)
    δ = (a.index == b.index)
    return δ + b*a
end

#σi*σj
function merge_idx_transitions(σ1::SymbolicUtils.Sym,σ2::SymbolicUtils.Sym)
    op = merge_idx_transitions(_to_qumulants(σ1), _to_qumulants(σ2))
    return _to_symbolic(op)
end
function merge_idx_transitions(σ1, σ2)
    i = σ1.index
    j = σ2.index
    δ = i==j
    σ3 = merge_transitions(σ1,σ2)
    if iszero(σ3)
        !(δ)*nip(σ1,σ2)
    else
        return δ*σ3[i] + !(δ)*nip(σ1,σ2)
    end
end

#nip*σ
function merge_nip_idx_transition(nip_::SymbolicUtils.Symbolic,σ::SymbolicUtils.Symbolic)
    p = merge_nip_idx_transition(_to_qumulants(nip_), _to_qumulants(σ))
    return _to_symbolic(p)
end
function merge_nip_idx_transition(nip_::OperatorTerm, σ::IndexedTransition)
    nip_σs = copy.(nip_.arguments)
    nip_idxs = find_index(nip_)
    if σ.index in nip_idxs
        it_σ = findfirst(x->isequal(x.index,σ.index), nip_σs)
        σ_new = merge_transitions(nip_σs[it_σ], σ)
        if iszero(σ_new)
            return 0
        else
            nip_σs[it_σ] = σ_new[σ.index]
            return nip(nip_σs)
        end
    else #σ.index not in nip_idxs
        all_summands = []
        δs = Number[]
        for itσ=1:length(nip_σs)
            σ_new = merge_transitions(nip_σs[itσ], σ)
            δ_ = (σ.index==nip_σs[itσ].index)
            push!(δs, !(δ_))
            if !iszero(σ_new)
                nip_σs_new = copy(nip_σs)
                nip_σs_new[itσ] = σ_new[σ.index]
                nip_σ_new = δ_*nip(nip_σs_new)
                push!(all_summands, nip_σ_new)
            end
        end
        push!(all_summands, *(δs...)*nip([nip_σs...,σ]))
        return sum(all_summands)
    end
end

#σ*nip
function merge_idx_transition_nip(σ::SymbolicUtils.Symbolic,nip_::SymbolicUtils.Symbolic)
    p = merge_idx_transition_nip( _to_qumulants(σ),_to_qumulants(nip_))
    return _to_symbolic(p)
end
function merge_idx_transition_nip(σ, nip_)
    nip_σs = copy(nip_.arguments)
    nip_idxs = find_index(nip_)
    if σ.index in nip_idxs
        it_σ = findfirst(x->isequal(x.index,σ.index), nip_σs)
        σ_new = merge_transitions(σ, nip_σs[it_σ])
        if iszero(σ_new)
            return 0
        else
            nip_σs[it_σ] = σ_new[σ.index]
            return nip(nip_σs)
        end
    else #σ.index not in nip_idxs
        all_summands = []
        δs = Number[]
        for itσ=1:length(nip_σs)
            σ_new = merge_transitions(σ,nip_σs[itσ])
            δ_ = (σ.index==nip_σs[itσ].index)
            push!(δs, !(δ_))
            if !iszero(σ_new)
                nip_σs_new = copy(nip_σs)
                nip_σs_new[itσ] = σ_new[σ.index]
                nip_σ_new = δ_*nip(nip_σs_new)
                push!(all_summands, nip_σ_new)
            end
        end
        push!(all_summands, *(δs...)*nip([σ,nip_σs...]))
        return sum(all_summands)
    end
end

# nip*nip TODO clean up
function merge_nips(nip1_::SymbolicUtils.Symbolic,nip2_::SymbolicUtils.Symbolic)
    p = merge_nips(_to_qumulants(nip1_), _to_qumulants(nip2_))
    return _to_symbolic(p)
end
function merge_nips(nip1_, nip2_)
    nip1_σs = copy(nip1_.arguments)
    nip2_σs = copy(nip2_.arguments)
    nip1_idxs = find_index(nip1_)
    nip2_idxs = find_index(nip2_)
    same_idxs = intersect(nip1_idxs, nip2_idxs)
    if isempty(same_idxs) #all indices unequal in both nips
        return merge_uneq_nips(nip1_σs, nip2_σs)
    else #identical indices in both nips, e.g. (σi⋅σj)*(σk⋅σj)
        same_idx_prods = []
        for idx in same_idxs
            nip1_σ_idx = nip1_σs[findfirst(x->isequal(x.index,idx),nip1_σs)]
            nip2_σ_idx = nip2_σs[findfirst(x->isequal(x.index,idx),nip2_σs)]
            δ_tmp = nip1_σ_idx.index==nip2_σ_idx.index
            σ_tmp = merge_transitions(nip1_σ_idx, nip2_σ_idx)
            iszero(σ_tmp) && return 0
            push!(same_idx_prods, σ_tmp[nip1_σ_idx.index])
        end
        same_idx_prods = sort(same_idx_prods, by=hash)
        nip1_σs_uneq = filter(x-> !(x.index in same_idxs), nip1_σs)
        nip2_σs_uneq = filter(x-> !(x.index in same_idxs), nip2_σs)

        if iszero(length(nip1_σs_uneq)*length(nip2_σs_uneq)) #all σs from one nip equal to the others
            nip_σs = sort(unique!(union(same_idx_prods, nip1_σs_uneq, nip2_σs_uneq)), by=hash)
            return nip(nip_σs)
        else
            return merge_uneq_nips(nip1_σs_uneq, nip2_σs_uneq;extra_σs=same_idx_prods)
        end
    end
end
function merge_uneq_nips(nip1_σs, nip2_σs; extra_σs=[]) #all indices unequal
    # extra_σs: needed for merge_nips (identical indices)
    len1 = length(nip1_σs)
    len2 = length(nip2_σs)
    len_prod = len1 + len2
    idx1 = Iterators.flatten(find_index.(nip1_σs))
    idx2 = Iterators.flatten(find_index.(nip2_σs))
    nip1_σs_ext = [nip1_σs;ones(Int,len2)] #extend array with "1"
    nip2_σs_ext = [nip2_σs;ones(Int,len1)]
    nip2_permu = unique(permutations(nip2_σs_ext))
    #
    # return [merge_eq_idx_trans(nip1_σs_ext[it], nip2_permu[1][it])[2] for it=1:len_prod]
    nip_sum_ls = []
    for it2=1:length(nip2_permu)
        check_zero = false
        δ_tmp_ls = []
        σ_tmp_ls = []
        for it1=1:len_prod
            isone(nip1_σs_ext[it1]) && isone(nip2_permu[it2][it1]) && continue
            if isone(nip1_σs_ext[it1])
                for ii_ in idx1
                    push!(δ_tmp_ls, nip2_permu[it2][it1].index != ii_)
                end
                push!(σ_tmp_ls, nip2_permu[it2][it1])
            elseif isone(nip2_permu[it2][it1])
                for ii_ in idx2
                    push!(δ_tmp_ls, nip1_σs_ext[it1].index != ii_)
                end
                push!(σ_tmp_ls, nip1_σs_ext[it1])
            else
                δ_tmp = nip1_σs_ext[it1].index==nip2_permu[it2][it1].index
                σ_tmp = merge_transitions(nip1_σs_ext[it1], nip2_permu[it2][it1])
                push!(δ_tmp_ls, δ_tmp)
                check_zero = iszero(σ_tmp)
                check_zero && break
                push!(σ_tmp_ls, σ_tmp[nip1_σs_ext[it1].index])
            end
            check_zero && break
        end
        check_zero && continue

        filter!(!isone, δ_tmp_ls)
        δ_tmp_ls = sort(δ_tmp_ls, by=hash)
        if isempty(δ_tmp_ls)
            δ_tmp_ls = [1]
        end
        filter!(!isone, σ_tmp_ls)
        sort!(σ_tmp_ls, by=hash)
        push!(nip_sum_ls, [δ_tmp_ls, σ_tmp_ls])
    end
    nip_sum_ls = unique!(nip_sum_ls)
    nip_sum = [prod(args[1])*nip(sort([args[2]...,extra_σs...], by=hash)) for args in nip_sum_ls]
    return sum(nip_sum)
end


# Rewriting σgg[i]
function rewrite_gs(σ::IndexedTransition)
    h = σ.hilbert
    aon = acts_on(σ)
    gs = ground_state(h,aon)
    i,j = σ.i, σ.j
    if i==j==gs
        idx = σ.index
        args = Any[1]
        for k in levels(h,aon)
            (k==i) || push!(args, -1*IndexedTransition(h, σ.name, k, k, aon, idx))
        end
        return +(args...)
    else
        return nothing
    end
end

### Parameters

struct IndexedParameter{T<:Number,S,I} <: SymbolicNumber{T}
    name::S
    index::I
    function IndexedParameter{T,S,I}(name::S,idx::I) where {T,S,I}
        param_idx = new(name,idx)
        if !haskey(IDX_TO_SYMS, param_idx)
            sym = SymbolicUtils.Sym{IndexedParameter}(gensym(:IndexedParameter))
            IDX_TO_SYMS[param_idx] = sym
            SYMS_TO_IDX[sym] = param_idx
        end
        return param_idx
    end
end
IndexedParameter{T}(name::S,idx::I) where {T,S,I} = IndexedParameter{T,S,I}(name,idx)
IndexedParameter(name, idx) = IndexedParameter{Number}(name, idx)
_to_symbolic(param_idx::IndexedParameter) = IDX_TO_SYMS[param_idx]
_to_qumulants(t::SymbolicUtils.Sym{T}) where T<:IndexedParameter = SYMS_TO_IDX[t]

Base.isequal(p::T, q::T) where T<:IndexedParameter = (p.name==q.name && isequal(p.index,q.index))
Base.hash(p::IndexedParameter{T}, h::UInt) where T = hash(p.name, hash(p.index, hash(T, h)))

# Methods
Base.conj(p::IndexedParameter{<:Real}) = p

find_index(s::SymbolicUtils.Symbolic) = _to_symbolic.(find_index(_to_qumulants(s)))
function find_index(t::OperatorTerm)
    idx = Index[]
    for arg in t.arguments
        append!(idx, find_index(arg))
    end
    unique!(idx)
    return idx
end
find_index(x) = Index[]
find_index(x::Union{IndexedTransition,IndexedCreate,IndexedDestroy}) = [x.index]
find_index(x::IndexedParameter) = x.index


### Symbolic Summation
Sum(ops::AbstractOperator, index::Index...) = OperatorTerm(Sum, [ops, index...])
Sum(x::SymbolicNumber, index::Index...) = NumberTerm(Sum, [x, index...])
# Sum(ops::Number, index::Index) = NumberTerm(Sum, [ops, index])
# Sum(ops::NumberTerm, index::Index) = NumberTerm(Sum, [ops, index])


Sum(op::SymbolicUtils.Symbolic{<:AbstractOperator}, index...) = SymbolicUtils.Term{AbstractOperator}(Sum, [op, index...])
Sum(op::Union{SymbolicUtils.Symbolic{<:Number},Number}, index...)= SymbolicUtils.Term{Number}(Sum, [op, index...])
# find_index(Sum_::OperatorTerm{typeof(Sum)}) = Sum_.arguments[2]

# swap index
swap_index(x, i::Union{Index,Int}, j::Union{Index,Int}) = x
swap_index(ex::Union{Index,Int}, i1::Union{Index,Int}, i2::Union{Index,Int}) = (isequal(ex,i1) ? i2 : ex)
swap_index(ex::IndexedTransition, i1::Union{Index,Int}, i2::Union{Index,Int}) = (isequal(ex.index,i1) ? IndexedTransition(ex.hilbert, ex.name, ex.i, ex.j, ex.aon, i2) : ex)
swap_index(ex::IndexedDestroy, i1::Union{Index,Int}, i2::Union{Index,Int}) = (isequal(ex.index,i1) ? IndexedDestroy(ex.hilbert, ex.name, ex.aon, i2) : ex)
swap_index(ex::IndexedCreate, i1::Union{Index,Int}, i2::Union{Index,Int}) = (isequal(ex.index,i1) ? IndexedCreate(ex.hilbert, ex.name, ex.aon, i2) : ex)
function swap_index(ex::IndexedParameter, i1::Union{Index,Int}, i2::Union{Index,Int})
    inds = Number[ex.index...]
    for idx in findall(isequal(i1), ex.index)
        inds[idx] = i2
    end
    return IndexedParameter(ex.name, inds)
end
swap_index(ex::OperatorTerm, i1::Union{Index,Int}, i2::Union{Index,Int}) = OperatorTerm(ex.f, [swap_index(arg, i1, i2) for arg in ex.arguments])
swap_index(ex::NumberTerm, i1::Union{Index,Int}, i2::Union{Index,Int}) = NumberTerm(ex.f, [swap_index(arg, i1, i2) for arg in ex.arguments])
swap_index(ex::SymbolicUtils.Symbolic, i::SymbolicUtils.Symbolic, j::SymbolicUtils.Symbolic) = _to_symbolic(swap_index(_to_qumulants(ex), _to_qumulants(i), _to_qumulants(j)))

function _multiply_idxs_borders(x, inds)
    args = Any[x]
    idx_ = _to_qumulants.(inds)
    for i in idx_
        push!(args, i.count+1)
    end
    return _to_symbolic(*(args...))
end

# TODO: write rule to extract non-indexed symbolic numbers out of sums
has_indexed(x) = false
has_indexed(s::SymbolicUtils.Symbolic) = has_indexed(_to_qumulants(s))
function has_indexed(t::NumberTerm)
    for arg in t.arguments
        has_indexed(arg) && return true
    end
    return false
end
has_indexed(::IndexedParameter) = true
has_indexed(::AbstractOperator) = true

sort_idx(idx) = SymbolicUtils.arguments(SymbolicUtils.sort_args(*, idx))

rewrite_nip_times(x) = x
function rewrite_nip_times(op::OperatorTerm)
    args = []
    for arg in op.arguments
        push!(args, rewrite_nip_times(arg))
    end
    return op.f(args...)
end
rewrite_nip_times(op::OperatorTerm{<:typeof(nip)}) = *(op.arguments...)
