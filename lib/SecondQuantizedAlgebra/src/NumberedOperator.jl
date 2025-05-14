#For representing in average terms
"""

    NumberedOperator <: QSym

Defines an operator, associated with a Number. Commutator-relations are calculated using these numbers, as a sort of a specific index-value.

Fields:
======

* op: An Operator, either a [`Transition`](@ref), a [`Destroy`](@ref) or a [`Create`](@ref) can be defined.
* numb: An Integer Number.

"""
struct NumberedOperator <:QSym
    op::IndexableOps
    numb::Int64
    function NumberedOperator(op::IndexableOps,numb::Int64)
        (numb <= 0) && error("can not create numbered-operator with negative or 0 number")
        return new(op,numb)
    end
end
function NumberedOperator(op,numb::Int64)
    (op isa SNuN) && return op
    if SymbolicUtils.iscall(op)
        f = SymbolicUtils.operation(op)
        args = [NumberedOperator(arg,numb) for arg in SymbolicUtils.arguments(op)]
        isempty(args) && return 0
        isequal(length(args),1) && return args[1]
        return f(args...)
    end
end

#define calculus for numbered operators -> break it down into QNumber multiplication

ismergeable(a::NumberedOperator,b::NumberedOperator) = isequal(a.numb,b.numb) ? ismergeable(a.op,b.op) : false

*(numOp::NumberedOperator, qmul::QMul) = inorder!(QMul(qmul.arg_c,vcat(numOp,qmul.args_nc)))
*(qmul::QMul, numOp::NumberedOperator) = inorder!(QMul(qmul.arg_c,vcat(qmul.args_nc,numOp)))

function *(numOp1::NumberedOperator,numOp2::NumberedOperator)
    if numOp1.op isa Create || numOp1.op isa Destroy || numOp2.op isa Create || numOp2.op isa Destroy
        return merge_commutators(1,[numOp1,numOp2])
    end
    return (numOp1.numb == numOp2.numb && isequal(acts_on(numOp1.op),acts_on(numOp2.op))) ? NumberedOperator(numOp1.op*numOp2.op,numOp1.numb) : QMul(1,sort([numOp1,numOp2],by=get_numbers))
end

hilbert(x::NumberedOperator) = hilbert(x.op)
Base.adjoint(x::NumberedOperator) = NumberedOperator(Base.adjoint(x.op),x.numb)
has_cluster(x::NumberedOperator) = has_cluster(x.op)
acts_on(x::NumberedOperator) = acts_on(x.op)

function Base.hash(a::NumberedOperator,h::UInt)
    return hash(NumberedOperator, hash(a.op, hash(a.numb, h)))
end
Base.isequal(x::NumberedOperator,y::NumberedOperator) = isequal(x.op,y.op) && isequal(x.numb,y.numb)
