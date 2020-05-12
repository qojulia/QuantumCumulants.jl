using Qumulants
using Test
using SymbolicUtils

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))

a_ = Destroy(hf,:a)
a = a_⊗one(ha)
σ_ = Transition(ha,:σ,:g,:e)
σ = one(hf)⊗σ_

ex = a*(a + a)
# @test simplify_operators(ex) == a_^2 ⊗ one(ha)

# JC
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

commutator(a,b) = simplify_operators(a*b - b*a)
da = commutator(1.0im*H,a)
@test da == -1.0im*ωc*a + (-1.0im*g)*σ
ds = commutator(1.0im*H,σ)
@test ds == (-1.0im*g)*a + (2.0im*g)*a_⊗Transition(ha,:σ,:e,:e) + (-1.0im*ωa)*σ



# a_sym = Qumulants._to_symbolic(a_)
# ad_sym = Qumulants._to_symbolic(a_')
# σ_sym = Qumulants._to_symbolic(σ_)
# σd_sym = Qumulants._to_symbolic(σ_')
# idc = Qumulants._to_symbolic(one(a_))
# ida = Qumulants._to_symbolic(one(σ_))
# zra = Qumulants._to_symbolic(zero(σ_))
#
# ex = 2*(σ'*a)*σ - σ*(2*(σ'*a))
# ex2 = simplify_operators(ex;fixpoint=false)
# ex3 = (((-2*a_) ⊗ Transition(ha,:σ,:g,:g))+((2*a_) ⊗ (σ_'*σ_)))
# ex4 = simplify_operators(ex3;rules=Qumulants.COMMUTATOR_RULES)
# ex5 = simplify_operators(simplify_operators(ex4;rules=Qumulants.TENSOR_RULES);rules=Qumulants.NC_TIMES_RULES)
#
# ex_sym = Qumulants._to_symbolic(ex4)
#
# ex1 = simplify_operators(ex.arguments[1])
# ex2 = simplify_operators(-1*ex.arguments[2])
#
# args = SymbolicUtils.Term[(-2*a_sym)⊗ida, (2*a_sym)⊗(SymbolicUtils.Term(*,[Qumulants._to_symbolic(Transition(ha,:σ,:e,:e))])), (-2*a_sym)⊗(-1*Qumulants._to_symbolic(Transition(ha,:σ,:e,:e)))]
