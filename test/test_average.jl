using QuantumCumulants
using SymbolicUtils
using Test

@testset "average" begin

    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf⊗ha

    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @test isequal(average(2a), 2*average(a))

    @test isequal(QuantumCumulants._conj(average(2im*a'*σ)), (-2im)*average(a*σ'))
    @test isequal(average(2*(a+a)), 2*(average(a) + average(a)))
    @test isequal(average(a^2), average(a*a))
    @test isequal(average(a^2, 1), average(a)^2)
    @test isequal(average((a'*a)^2), average(a'*a*a'*a))
    @test isequal(
        average(a'*a^2, 2),
        average(a')*average(a^2) + -2*average(a')*average(a)^2 + 2*average(a)*average(a'*a),
    )

    @test isequal(simplify(average(σ)+average(σ)), average(2σ))
    @test iszero(simplify(average(σ)-average(σ)))

    ωc, ωa = cnumbers("ω_c ω_a")
    @test isequal(average(ωc), ωc)
    @test isequal(average(ωc*a), ωc*average(a))
    @test isequal(average(ωc*(a+a')), ωc*average(a) + ωc*average(a'))

    n = average(a'*a)
    @test isequal(cumulant_expansion(n, 2), n)
    @test isequal(cumulant_expansion(n, 1), average(a')*average(a))

    @test iszero(average(a*σ)*average(a) - average(a)*average(a*σ))

    # Test cumulants
    hs = FockSpace[]
    for i = 1:4
        push!(hs, FockSpace(Symbol(:fock, i)))
    end
    h = ⊗(hs...)

    a = Destroy(h, :a, 1)
    b = Destroy(h, :b, 2)
    c = Destroy(h, :c, 3)
    d = Destroy(h, :d, 4)

    @test isequal(cumulant(a*b), simplify(average(a*b)+-1*average(a)*average(b)))
    @test isequal(cumulant(a*b, 1), average(a*b))
    @test iszero(
        simplify(
            cumulant(a*b*c) - (
                average(a*b*c) + 2*average(a)*average(b)*average(c) -
                average(a)*average(b*c) - average(b)*average(a*c) - average(c)*average(a*b)
            ),
        ),
    )
    @test isequal(
        simplify(average(a*b*c*d) - cumulant(a*b*c*d)),
        cumulant_expansion(average(a*b*c*d), 3),
    )

@test get_order(a*b + c) == 2

end # testset
