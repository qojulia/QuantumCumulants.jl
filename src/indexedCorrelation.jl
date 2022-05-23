include("indexedScale.jl")

_new_operator(op::IndexedOperator,h,aon=acts_on(op)) = IndexedOperator(Transition(h,op.op.name,op.op.i,op.op.j,aon;op.op.metadata),op.ind)
_new_operator(nOp::NumberedOperator,h,aon=acts_on(nOp)) = NumberedOperator(Transition(h,nOp.op.name,nOp.op.i,nOp.op.j,aon;nOp.op.metadata),nOp.numb)