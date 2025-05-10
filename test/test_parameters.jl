using QuantumCumulants
using SymbolicUtils
using Test

@testset "parameters" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf⊗ha

    a = Destroy(h, :a)

    @cnumbers p q
    # @test p*q isa NumberTerm
    @test p*q*a isa QuantumCumulants.QTerm

    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e)

    # JC
    @cnumbers ωc ωa g
    H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

    da = commutator(1im*H, a)
    @test iszero(simplify(da - ((-1im*g)*σ + (-1im*ωc)*a)))
    ds = commutator(1im*H, σ)
    @test iszero(simplify(ds - ((-1im*g)*a + (-1im*ωa)*σ + (2im*g)*a*σee)))
    dn = commutator(1im*H, a'*a)
    @test iszero(simplify(dn - ((-1im*g)*a'*σ + (1im*g)*a*σ')))


end # testset
