using Qumulants
using SymbolicUtils
SymbolicUtils.show_simplified[]=false
using Test

# Test
hf = FockSpace(:c)

a = Destroy(hf,:a)

h2 = hf⊗hf
a_embed = embed(h2,a,2)
a_sym = Qumulants._to_symbolic(a_embed)
@test Qumulants._to_qumulants(a_sym)==a_embed

A(i) = embed(h2,a,i)

@test A(1)*A(2)'==simplify_operators(A(1)*A(2)')

# Composite system
commutator(a,b) = simplify_operators(a*b - b*a)

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = embed(h,Destroy(hf,:a),1)
σ = embed(h,Transition(ha,:σ,:g,:e),2)

@test commutator(a'*a,a) == -1*a
@test commutator(σ'*σ,σ) == -1*σ

@test commutator(1*a'*a,a) == -1*a
@test commutator(1*σ'*σ,σ) == -1*σ

tmp = commutator(1*σ'*σ,σ)
ex = (-1 + embed(h,Transition(ha,:σ,:e,:e),2))*σ

# JC Model
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1.0im*H,a)
@test da == -1.0im*ωc*a + (-1.0im*g)*σ
ds = commutator(1.0im*H,σ)
@test ds == ((0.0-1.0im)*g)*a + ((0.0-1.0im)*ωa)*σ  + (2.0im*g)*a*embed(h,Transition(ha,:σ,:e,:e),2)
