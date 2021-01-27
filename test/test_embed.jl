using Qumulants
using Test

@testset "embed" begin

hf = FockSpace(:c)

a = Destroy(hf,:a)

h2 = hf⊗hf
a_embed = embed(h2,a,2)
@test_throws ErrorException Destroy(h2,:a)
@test isequal(a_embed , Destroy(h2,:a,2))

a1 = Destroy(h2,:a,1)
a2 = Destroy(h2,:a,2)
@test isequal(a1*a2',qsimplify(a1*a2'))

# Composite system
commutator(a,b) = qsimplify(a*b - b*a)

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

@test isequal(embed(h,Destroy(hf,:a),1),a)
@test isequal(embed(h,Transition(ha,:σ,:g,:e),2),σ)

@test_throws AssertionError Destroy(h,:a,2)
@test_throws AssertionError Create(h,:a,2)
@test_throws AssertionError Transition(h,:σ,:g,:e,1)

@test isequal(commutator(a'*a,a) , -1*a)
@test isequal(commutator(σ'*σ,σ) , -1*σ)

@test isequal(commutator(1*a'*a,a) , -1*a)
@test isequal(commutator(1*σ'*σ,σ) , -1*σ)

tmp = commutator(1*σ'*σ,σ)
ex = (-1 + embed(h,Transition(ha,:σ,:e,:e),2))*σ

# JC Model
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1.0im*H,a)
@test iszero(qsimplify(da - (-1.0im*ωc*a + (-1.0im*g)*σ)))
ds = commutator(1.0im*H,σ)
@test iszero(qsimplify(ds - (((0.0-1.0im)*g)*a + ((0.0-1.0im)*ωa)*σ  + (2.0im*g)*a*Transition(h,:σ,:e,:e))))

end # testset
