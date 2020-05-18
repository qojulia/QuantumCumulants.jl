abstract type SymbolicNumber <: Number end

struct NumberTerm{T<:Number,F,ARGS} <: SymbolicNumber
    f::F
    arguments::ARGS
end
function NumberTerm(f,args)
    T = promote_type(typeof.(args)...)
    return NumberTerm{T}(f,args)
end
Base.:(==)(t1::NumberTerm,t2::NumberTerm) = (t1.f===t2.f && t1.arguments==t2.arguments)

Base.zero(::Type{<:SymbolicNumber}) = 0
Base.one(::Type{<:SymbolicNumber}) = 1
Base.conj(t::NumberTerm) = NumberTerm(t.f, conj.(t.arguments))

struct Parameter{T<:Number} <: SymbolicNumber
    name::Symbol
end
Parameter(name::Symbol) = Parameter{Number}(name)

Base.conj(p::Parameter{<:Real}) = p
Base.conj(p::Parameter) = NumberTerm(conj, [p])

for f = [:+,:-,:*,:/,:^]
    @eval Base.$f(a::SymbolicNumber,b::Number) = NumberTerm($f, [a,b])
    @eval Base.$f(a::Number,b::SymbolicNumber) = NumberTerm($f, [a,b])
    @eval Base.$f(a::SymbolicNumber,b::SymbolicNumber) = NumberTerm($f, [a,b])
end
for f = [:cos,:sin,:tan,:sqrt]
    @eval Base.$f(a::SymbolicNumber) = NumberTerm($f, [a])
end

# Conversion to SymbolicUtils
_to_symbolic(p::Parameter{T}) where T = SymbolicUtils.Sym{T}(p.name)
_to_symbolic(n::NumberTerm{T}) where T = SymbolicUtils.Term{T}(n.f, _to_symbolic.(n.arguments))
function _to_qumulants(s::SymbolicUtils.Sym{T}) where T<:Number
    # if haskey(SYMS_TO_AVERAGES, s)
    #     return SYMS_TO_AVERAGES[s]
    # else
    return Parameter{T}(s.name)
    # end
end
_to_qumulants(t::SymbolicUtils.Term{T}) where T<:Number = NumberTerm{T}(t.f, _to_qumulants.(t.arguments))


# macro parameters(xs...)
#     ex = Expr(:block)
#     pnames = Symbol[]
#     for x in xs
#         push!(pnames, x)
#         ex_ = :( Parameter($(Meta.quot(x))) )
#         push!(ex.args, :( $x = $ex_ ) )
#     end
#     # push!(ex.args, Expr(:tuple, pnames...))
#     return ex
# end


# macro symbolic(ex)
#     @assert ex.head === :(=)
#     lhs = ex.args[1]
#     rhs = ex.args[2]
#     @assert rhs.args[1]∈INTERFACE_OPERATORS
#     if rhs.args[1]==:destroy
#         hname = Expr(:quote, rhs.args[2])
#         h = :( FockSpace($hname) )
#         return :($(esc(lhs)) = Destroy($h, $(Expr(:quote, lhs))) )
#     elseif rhs.args[1]==:transition
#     end
# end
