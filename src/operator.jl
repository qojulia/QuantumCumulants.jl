import SymbolicUtils
# Abstract types
abstract type AbstractOperator end
abstract type BasicOperator <: AbstractOperator end

isoperator(x) = false
isoperator(x::Union{T,SymbolicUtils.Symbolic{T}}) where T<:AbstractOperator = true

# Keep track of operators when converting to SymbolicUtils
const OPERATORS_TO_SYMS = Dict{BasicOperator,SymbolicUtils.Sym}()
const SYMS_TO_OPERATORS = Dict{SymbolicUtils.Sym,BasicOperator}()

# Operator expressions
struct OperatorTerm{F,ARGS} <: AbstractOperator
    f::F
    arguments::ARGS
end
Base.:(==)(t1::OperatorTerm,t2::OperatorTerm) = (t1.f===t2.f && t1.arguments==t2.arguments)

for f = [:+,:-,:*]
    @eval Base.$f(a::AbstractOperator,b::AbstractOperator) = (check_hilbert(a,b); OperatorTerm($f, [a,b]))
end
Base.:*(a::AbstractOperator,b::Number) = OperatorTerm(*, [a,b])
Base.:*(a::Number,b::AbstractOperator) = OperatorTerm(*, [a,b])
Base.:^(a::AbstractOperator,b) = OperatorTerm(^, [a,b])

# Variadic methods
Base.:+(x::AbstractOperator) = x
Base.:-(x::AbstractOperator) = -1*x
Base.:+(x::AbstractOperator, w::AbstractOperator...) = (rec_check_hilbert(x,w); OperatorTerm(+, [x;w...]))

Base.:*(x::AbstractOperator) = x
Base.:*(x::AbstractOperator, w...) = (rec_check_hilbert(x,w); OperatorTerm(*, [x;w...]))
Base.:*(x, y::AbstractOperator, w...) = (rec_check_hilbert(x,y,w); OperatorTerm(*, [x;y;w...]))
Base.:*(x::AbstractOperator, y::AbstractOperator, w...) = (rec_check_hilbert(x,y,w); OperatorTerm(*, [x;y;w...]))

Base.adjoint(t::OperatorTerm) = OperatorTerm(t.f, adjoint.(t.arguments))
Base.adjoint(t::OperatorTerm{<:typeof(*)}) = OperatorTerm(t.f, reverse(adjoint.(t.arguments)))

check_hilbert(a,b) = nothing # TODO
rec_check_hilbert(args...) = nothing

# Tensor product
⊗(a::AbstractOperator,b::AbstractOperator) = OperatorTerm(⊗, [a,b])
⊗(a::SymbolicUtils.Symbolic{<:AbstractOperator},b::SymbolicUtils.Symbolic{<:AbstractOperator}) = SymbolicUtils.Term(⊗, [a,b])

# Variadic ⊗
⊗(a::AbstractOperator) = a
⊗(a::AbstractOperator,b::AbstractOperator...) = OperatorTerm(⊗, [a;b...])
⊗(a::SymbolicUtils.Symbolic{<:AbstractOperator},b::SymbolicUtils.Symbolic{<:AbstractOperator}...) = SymbolicUtils.Term(⊗, [a;b...])

# Function hierarchy * < ⊗ < +
import Base: *, +
fs = [:⊗, :*]
for f in fs
    @eval $(f)(a::AbstractOperator, b::OperatorTerm{<:typeof(+)}) = OperatorTerm(+, [$(f)(a,b_) for b_ in b.arguments])
    @eval $(f)(a::OperatorTerm{<:typeof(+)}, b::AbstractOperator) = OperatorTerm(+, [$(f)(a_,b) for a_ in a.arguments])
    @eval $(f)(a::OperatorTerm{<:typeof(+)}, b::OperatorTerm{<:typeof(+)}) = OperatorTerm(+, [$(f)(a_,b_) for a_ in a.arguments, b_ in b.arguments][:])
end
*(a::Number, b::OperatorTerm{<:typeof(+)}) = OperatorTerm(+, [*(a,b_) for b_ in b.arguments])
function *(a::OperatorTerm{<:typeof(⊗)}, b::OperatorTerm{<:typeof(⊗)})
    @assert length(a.arguments)==length(b.arguments)
    OperatorTerm(⊗, [*(a_,b_) for (a_,b_) in zip(a.arguments,b.arguments)])
end
*(a::Number, b::OperatorTerm{<:typeof(⊗)}) = OperatorTerm(⊗, [*(a,b.arguments[1]); b.arguments[2:end]])
Base.:^(a::OperatorTerm{<:typeof(⊗)}, b) = OperatorTerm(⊗, [^(a_,b) for a_ in a.arguments])

# General basic operator types
struct Identity{H,S} <: BasicOperator
    hilbert::H
    name::S
    function Identity{H,S}(hilbert::H,name::S) where {H,S}
        op = new(hilbert,name)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = SymbolicUtils.Sym{Identity}(gensym(:Identity))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Identity(hilbert::H,name::S) where {H,S} = Identity{H,S}(hilbert,name)
Identity(hilbert::HilbertSpace) = Identity(hilbert, :𝟙)
Base.one(hilbert::HilbertSpace) = Identity(hilbert)
Base.one(a::BasicOperator) = one(a.hilbert)
Base.isone(::AbstractOperator) = false
Base.isone(::Identity) = true
Base.adjoint(x::Identity) = x
isidentity(x) = false
isidentity(x::Union{T,SymbolicUtils.Sym{T}}) where T<:Identity = true

struct Zero{H,S} <: BasicOperator
    hilbert::H
    name::S
    function Zero{H,S}(hilbert::H,name::S) where {H,S}
        op = new(hilbert,name)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = SymbolicUtils.Sym{Zero}(gensym(:Zero))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Zero(hilbert::H,name::S) where {H,S} = Zero{H,S}(hilbert,name)
Zero(hilbert::HilbertSpace) = Zero(hilbert, 0)
Base.zero(hilbert::HilbertSpace) = Zero(hilbert)
Base.zero(a::BasicOperator) = zero(a.hilbert)
Base.iszero(::AbstractOperator) = false
Base.iszero(::Zero) = true
Base.adjoint(x::Zero) = x
