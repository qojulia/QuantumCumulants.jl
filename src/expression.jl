import Base: ==, adjoint, *, +, -, ^

"""
    Expression{fType,argType} <: AbstractOperator

Type defining an expression which connects operators contained in its arguments
`args` by the function contained in `f`.

# Arguments
*`f`: The function which is applied to `args` `(*,+,⊗)`.
*`args`: Vector of expressions to which the derivatives are equal.

# Fields
*same as Arguments.
"""
mutable struct Expression{fType,argType} <: AbstractOperator
    f::fType
    args::argType
end
Expression(f,args...) = Expression(f,args)
# ==(a::Expression,b::Expression) = Expression(==,[a,b])
==(a::Expression,b::Expression) = (a.f==b.f && a.args==b.args)
Base.isapprox(a::Expression,b::Expression) = (a.f==b.f && isapprox(a.args,b.args))
# TODO: Add general flag to Expression determining whether there are indexed objects in the expression

const Add{argType} = Expression{typeof(+),argType}
const Prod{argType} = Expression{typeof(*),argType}

"""
    ⊗(a::AbstractOperator,b::AbstractOperator)

Compute the tensor product of two operators or expressions `a` and `b`.
"""
function ⊗ end
const TensorProd{argType} = Expression{typeof(⊗),argType}

function ==(a::Add,b::Add)
    length(a.args) == length(b.args) || return false
    check = true
    for i=1:length(a.args)
        check_ = false
        for j=1:length(b.args)
            check_ = (a.args[i] == b.args[j])
            check_ && break
        end
        check = check_
        check || break
    end
    return check
end

copy(a::Expression) = Expression(a.f, copy.(a.args))
adjoint(a::Expression) = Expression(a.f, adjoint.(a.args))
function adjoint(a::Prod)
    c = 1
    args = []
    for a1=a.args
        if isa(a1,Number)
            c *= conj(a1)
        else
            push!(args,adjoint(a1))
        end
    end
    reverse!(args)
    if isone(c)
        args_ = args
    else
        args_ = [c;args]
    end
    return Expression(*,args_)
end
Base.iszero(a::Union{Prod,TensorProd}) = any(iszero.(a.args))
Base.iszero(a::Add) = all(iszero.(a.args))
Base.one(a::Prod) = one(a.args[end])
Base.zero(a::Prod) = zero(a.args[end])
Base.one(a::Add) = one(a.args[1])
Base.zero(a::Add) = zero(a.args[1])
Base.one(a::TensorProd) = ⊗([one(a1) for a1=a.args]...)
Base.zero(a::TensorProd) = ⊗([zero(a1) for a1=a.args]...)

###########
# Algebra #
###########

# * Numbers
*(x::Number,a::BasicOperator) = Expression(*,[x,a])
*(a::BasicOperator,x::Number) = Expression(*,[x,a])

*(a::Number,b::Prod) = Expression(*,[a;b.args])
*(a::Prod,b::Number) = Expression(*,[a.args;b])

# * Operators
*(a::BasicOperator,b::BasicOperator) = Expression(*,[a,b])
*(a::BasicOperator,b::Prod) = Expression(*,[a;b.args])
*(a::Prod,b::BasicOperator) = Expression(*,[a.args;b])
*(a::Prod,b::Prod) = Expression(*,[a.args;b.args])

*(::Identity,a::Prod) = a
*(a::Prod,::Identity) = a
*(a::Zero,::Prod) = a
*(::Prod,a::Zero) = a


# ⊗
⊗(a::AbstractOperator,b::AbstractOperator) = Expression(⊗,[a,b])
⊗(a::AbstractOperator,b::TensorProd) = Expression(⊗,[a;b.args])
⊗(a::TensorProd,b::AbstractOperator) = Expression(⊗,[a.args;b])
⊗(a::TensorProd,b::TensorProd) = Expression(⊗,[a.args;b.args])
⊗(args...) = reduce(⊗, args)


*(x::Number,a::TensorProd) = Expression(⊗,[x*a.args[1]; a.args[2:end]])
*(a::TensorProd,x::Number) = x*a
function *(a::TensorProd,b::TensorProd)
    length(a.args) == length(b.args) || error("Tensor products must have the same number of arguments in order to be multiplied!")
    # return Expression(⊗,[a.args[i]*b.args[i] for i=1:length(a.args)])
    return ⊗([a.args[i]*b.args[i] for i=1:length(a.args)]...)
end

# +
+(a::AbstractOperator,b::AbstractOperator) = Expression(+,[a,b])
+(a::AbstractOperator,b::Add) = Expression(+,[a;b.args])
+(a::Add,b::AbstractOperator) = Expression(+,[a.args;b])
+(a::Add,b::Add) = Expression(+,[a.args;b.args])
function +(a::TensorProd,b::TensorProd)
    length(a.args) == length(b.args) || error("Tensor products must have the same number of arguments in order to be added!")
    return Expression(+,[a,b])
end
+(a::TensorProd,b::BasicOperator) = throw(MethodError(+,a,b))
+(a::BasicOperator,b::TensorProd) = +(b,a)
+(a::Prod,b::TensorProd) = throw(MethodError(+,a,b))
+(a::TensorProd,b::Prod) = +(b,a)
+(a::Prod,::Zero) = a
+(::Zero,a::Prod) = a

*(x::Number,a::Add) = Expression(+,[x*a1 for a1=a.args])
*(a::Add,x::Number) = x*a
*(a::AbstractOperator,b::Add) = Expression(+,[a*b1 for b1=b.args])
*(a::Add,b::AbstractOperator) = Expression(+,[a1*b for a1=a.args])
function *(a::Add,b::Add)
    args = []#typejoin(eltype(a.args),eltype(b.args))[]
    for a1=a.args, b1=b.args
        push!(args, a1*b1)
    end
    return sum(args)#Expression(+,args)
end

-(a::AbstractOperator) = -1*a
-(a::AbstractOperator,b::AbstractOperator) = a + (-b)

# ^
^(a::AbstractOperator,n::Int) = prod([a for i=1:n])

⊗(a::AbstractOperator,b::Add) = Expression(+,[a⊗b1 for b1=b.args])
⊗(a::Add,b::AbstractOperator) = Expression(+,[a1⊗b for a1=a.args])
⊗(a::TensorProd,b::Add) = Expression(+,[a⊗b1 for b1=b.args])
⊗(a::Add,b::TensorProd) = Expression(+,[a1⊗b for a1=a.args])
function ⊗(a::Add,b::Add)
    args = typejoin(eltype(a.args),eltype(b.args))[]
    for a1=a.args, b1=b.args
        push!(args, a1⊗b1)
    end
    return Expression(+,args)
end


function embed(a::BasicOperator,i::Int,n::Int)
    @assert 0 < i <= n
    ops = typejoin(typeof(a),Identity)[Identity() for j=1:n]
    ops[i] = a
    return ⊗(ops...)
end
⊗(a::AbstractOperator) = a
function embed(a::Vector{<:AbstractOperator},i::Vector{<:Int},n::Int)
    ops = typejoin(eltype(a),Identity)[Identity() for j=1:n]
    ops[i] = a
    return ⊗(ops...)
end

ishermitian(a::TensorProd) = all(ishermitian.(a.args))


# Sometimes it is necessary to avoid simplification
dont_simplify(ex::AbstractOperator) = Expression(dont_simplify, [ex])
dont_simplify(x::Number) = x
const DontSimplify{ARGS} = Expression{typeof(dont_simplify),ARGS}
dont_simplify(ex::DontSimplify) = ex

remove_dontsimplify(x::Number) = x
remove_dontsimplify(a::DontSimplify) = a.args[1]
remove_dontsimplify(a::AbstractOperator) = a
remove_dontsimplify(ex::Expression) = ex.f(remove_dontsimplify.(ex.args)...)

==(a::DontSimplify,b::DontSimplify) = (a.args==b.args)
==(a::AbstractOperator,b::DontSimplify) = (a==b.args[1])
==(a::DontSimplify,b::AbstractOperator) = (a.args[1]==b)
==(a::Expression,b::DontSimplify) = (a==b.args[1])
==(a::DontSimplify,b::Expression) = (a.args[1]==b)

Base.one(ex::DontSimplify) = one(ex.args[1])
Base.zero(ex::DontSimplify) = zero(ex.args[1])
Base.iszero(ex::DontSimplify) = iszero(ex.args[1])
Base.isone(ex::DontSimplify) = isone(ex.args[1])

*(x::Number,ex::DontSimplify) = dont_simplify(x*ex.args[1])
*(ex::DontSimplify,x::Number) = x*ex
*(ex::DontSimplify,a::BasicOperator) = Expression(*,[ex,a])
*(a::BasicOperator,ex::DontSimplify) = Expression(*,[a,ex])
*(a::DontSimplify,b::DontSimplify) = dont_simplify(a.args[1]*b.args[1])
*(ex::DontSimplify,a::Prod) = Expression(*,[ex,a])
*(a::Prod,ex::DontSimplify) = Expression(*,[a,ex])
