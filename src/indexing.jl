const IDX_TO_SYMS = Dict{SymbolicNumber,SymbolicUtils.Sym}()
const SYMS_TO_IDX = Dict{SymbolicUtils.Sym,SymbolicNumber}()

### Indices

struct Index{T<:Int,S,L,U} <: SymbolicNumber{T}
    name::S
    lower::L
    upper::U
    function Index{T,S,L,U}(name::S,lower::L,upper::U) where {T<:Int,S,L,U}
        idx = new(name,lower,upper)
        if !haskey(IDX_TO_SYMS, idx)
            sym = SymbolicUtils.Sym{Index}(gensym(:Index))
            IDX_TO_SYMS[idx] = sym
            SYMS_TO_IDX[sym] = idx
        end
        return idx
    end
end
Index{T}(name::S,lower::L,upper::U) where {T,S,L,U} = Index{T,S,L,U}(name,lower,upper)
Index(name,lower,upper) = Index{Int}(name,lower,upper)
_to_symbolic(idx::Index) = IDX_TO_SYMS[idx]
_to_qumulants(t::SymbolicUtils.Sym{T}) where T<:Index = SYMS_TO_IDX[t]

Base.hash(i::Index, h::UInt) = hash(i.upper, hash(i.lower, hash(i.name, h)))
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
nip(args::IndexedTransition...) = OperatorTerm(nip, args)

for T in [:Destroy,:Create,:Transition]
    Name = Symbol(:Indexed,T)
    @eval get_index(a::$Name) = a.index
end

### Simplification functions
function commute_bosonic_idx(a::SymbolicUtils.Symbolic,b::SymbolicUtils.Symbolic)
    return _to_symbolic(commute_bosonic_idx(_to_qumulants(a), _to_qumulants(b)))
end
function commute_bosonic_idx(a,b)
    δ = (a.index == b.index)
    return δ + b*a
end


### Parameters

struct IndexedParameter{T<:Number,S,I} <: SymbolicNumber{T}
    name::S
    idx::I
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

Base.isequal(p::T, q::T) where T<:IndexedParameter = (p.name==q.name && isequal(p.idx,q.idx))
Base.hash(p::IndexedParameter{T}, h::UInt) where T = hash(p.name, hash(p.idx, hash(T, h)))

# Methods
Base.conj(p::IndexedParameter{<:Real}) = p
