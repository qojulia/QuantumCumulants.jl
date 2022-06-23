using QuantumCumulants
using QuantumOpticsBase
using Test
using Random; Random.seed!(0)

@testset "numeric-conversion" begin

# Test fock basis conversion
hfock = FockSpace(:fock)
@qnumbers a::Destroy(hfock)
bfock = FockBasis(7)
@test to_numeric(a, bfock) == destroy(bfock)
@test to_numeric(a', bfock) == create(bfock)

# NLevelSpace conversion
hnlevel = NLevelSpace(:nlevel, 3)
σ(i,j) = Transition(hnlevel, :σ, i, j)
bnlevel = NLevelBasis(3)
for i=1:3, j=1:3
    op = σ(i,j)
    @test to_numeric(op, bnlevel) == transition(bnlevel, i, j)
end

# with symbolic levels
levels = (:g, :e, :a)
hnlevel_sym = NLevelSpace(:nlevel_sym, levels)
σ_sym(i,j) = Transition(hnlevel_sym, :σ, i, j)
@test_throws ArgumentError to_numeric(σ_sym(:e,:g), bnlevel)
level_map = Dict((levels .=> (1,2,3))...)
for i=1:3, j=1:3
    lvl1 = levels[i]
    lvl2 = levels[j]
    op = σ_sym(lvl1, lvl2)
    @test to_numeric(op, bnlevel; level_map=level_map) == transition(bnlevel, i, j)
end


# On composite bases
hprod = hfock ⊗ hnlevel
a = Destroy(hprod, :a)
σprod(i,j) = Transition(hprod, :σ, i, j)
bprod = bfock ⊗ bnlevel
for i=1:3, j=1:3
    op1 = a*σprod(i,j)
    op2 = a'*σprod(i,j)
    @test to_numeric(op1, bprod) == destroy(bfock) ⊗ transition(bnlevel, i, j)
    @test to_numeric(op2, bprod) == create(bfock) ⊗ transition(bnlevel, i, j)
end

@test to_numeric(a'*a, bprod) ≈ number(bfock) ⊗ one(bnlevel)

# Composite basis with symbolic levels
σsym_prod(i,j) = Transition(hfock ⊗ hnlevel_sym, :σ, i, j)
a = Destroy(hfock ⊗ hnlevel_sym, :a)
@test_throws ArgumentError to_numeric(a*σsym_prod(:e,:g), bprod)
for i=1:3, j=1:3
    op1 = a*σsym_prod(levels[i],levels[j])
    op2 = a'*σsym_prod(levels[i],levels[j])
    @test to_numeric(op1, bprod; level_map=level_map) == destroy(bfock) ⊗ transition(bnlevel, i, j)
    @test to_numeric(op2, bprod; level_map=level_map) == create(bfock) ⊗ transition(bnlevel, i, j)
end

# Numeric average values
α = 0.1 + 0.2im
ψ = coherentstate(bfock, α)
a = Destroy(hfock, :a)
@test numeric_average(a, ψ) ≈ numeric_average(average(a), ψ) ≈ α
@test numeric_average(a'*a, ψ) ≈ numeric_average(average(a'*a), ψ) ≈ abs2(α)

ψprod = ψ ⊗ nlevelstate(bnlevel, 1)
@test_throws ArgumentError numeric_average(σsym_prod(:e,:g), ψprod)
idfock = one(bfock)
for i=1:3, j=1:3
    op = σprod(i, j)
    op_sym = σsym_prod(levels[i],levels[j])
    op_num = idfock ⊗ transition(bnlevel, i, j)
    @test numeric_average(op, ψprod) ≈ expect(op_num, ψprod)
    @test numeric_average(op_sym, ψprod; level_map=level_map) ≈ expect(op_num, ψprod)
end


# Initial values in actual equations
levels = (:g,:e)
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, levels)
a = Destroy(h,:a)
s(i,j) = Transition(h, :σ, i, j)

@cnumbers Δ g κ γ η

H = Δ*a'*a + g*(a'*s(:g,:e) + a*s(:e,:g)) + η*(a + a')
ops = [a,s(:g,:e),a'*a,s(:e,:e),a'*s(:g,:e)]
eqs = meanfield(ops,H,[a];rates=[κ],order=2)

bcav = FockBasis(10)
batom = NLevelBasis(2)
b = bcav ⊗ batom
ψ0 = randstate(b)

level_map = Dict((levels .=> [1,2])...)
u0 = initial_values(eqs, ψ0; level_map=level_map)

@test u0[1] ≈ expect(destroy(bcav) ⊗ one(batom), ψ0)
@test u0[2] ≈ expect(one(bcav) ⊗ transition(batom, 1, 2), ψ0)
@test u0[3] ≈ expect(number(bcav) ⊗ one(batom), ψ0)
@test u0[4] ≈ expect(one(bcav) ⊗ transition(batom, 2, 2), ψ0)
@test u0[5] ≈ expect(create(bcav) ⊗ transition(batom, 1, 2), ψ0)


end # testset