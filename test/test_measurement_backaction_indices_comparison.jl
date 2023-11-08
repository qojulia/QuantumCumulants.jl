using QuantumCumulants
using Test
using SymbolicUtils
using Symbolics

@testset "test_measurement_backaction_indices_comparison" begin

@cnumbers N ωa γ η χ ωc κ g ξ ωl
@syms t::Real
@register pulse(t)

hc = FockSpace(:resonator)
ha = NLevelSpace(:atom,2)
h = hc ⊗ ha

j = Index(h,:j,N,ha)
k = Index(h,:k,N,ha)

@qnumbers a::Destroy(h,1)
σ(α,β,k) = IndexedOperator(Transition(h,:σ,α,β,2), k)

H = ωc * a' * a + ωa * Σ(σ(2,2,j),j)+g*a'*Σ(σ(1,2,j),j)+g*a*Σ(σ(2,1,j),j)
J = [a*exp(1.0im*ωl*t),σ(1,2,j),σ(2,1,j),σ(2,2,j)]
rates = [κ,γ,η*pulse(t),2*χ]
efficiencies = [ξ,0,0,0]
ops = [a,a'*a,σ(2,2,k),σ(1,2,k)]

eqs = meanfield(ops,H,J; rates = rates, efficiencies = efficiencies, order = 2)
eqs_c = complete(eqs)
scaled_eqs = scale(eqs_c)


operators = eqs_c.operators
full_stoch_eqs = meanfield(operators,H,J; rates = rates, efficiencies = efficiencies, order = 4)
full_det_eqs = meanfield(operators,H,J; rates = rates, order = 4)
stoch_eqs = meanfield(ops,H,J; rates = rates, efficiencies = efficiencies, order = 2)
det_eqs = meanfield(ops,H,J; rates = rates, order = 2)

stoch_eqs_complete = indexed_complete(stoch_eqs)
det_eqs_complete = indexed_complete(det_eqs)

for (eq_stoch, eq_det) in zip(full_stoch_eqs.equations, full_det_eqs.equations)
    @test isequal(simplify(eq_stoch.lhs - eq_det.lhs),0)
    @test isequal(simplify(eq_stoch.rhs - eq_det.rhs),0)
end

for (eq_stoch, eq_det) in zip(stoch_eqs.equations, det_eqs.equations)
    @test isequal(simplify(eq_stoch.lhs - eq_det.lhs),0)
    @test isequal(simplify(eq_stoch.rhs - eq_det.rhs),0)
end

test_mf = QuantumCumulants.IndexedMeanfieldEquations(full_stoch_eqs.noise_equations, full_stoch_eqs.operator_equations, full_stoch_eqs.states, full_stoch_eqs.operators, full_stoch_eqs.hamiltonian,full_stoch_eqs.jumps, full_stoch_eqs.jumps_dagger, full_stoch_eqs.rates, full_stoch_eqs.iv, full_stoch_eqs.varmap, 4)
test_mf = QuantumCumulants.cumulant_expansion(test_mf,2)

for (lhs, rhs) in zip(stoch_eqs_complete.noise_equations, test_mf.equations)
    @test isequal(simplify(lhs.lhs - rhs.lhs),0)
    @test isequal(simplify(lhs.rhs - rhs.rhs),0)
end


stoch_eqs_scaled = scale(stoch_eqs)
det_eqs_scaled = scale(det_eqs)

for (lhs, rhs) in zip(stoch_eqs_scaled.equations, det_eqs_scaled.equations)
    @test isequal(simplify(lhs.lhs - rhs.lhs),0)
    @test isequal(simplify(lhs.rhs - rhs.rhs),0)
end

test_mf = QuantumCumulants.IndexedMeanfieldEquations(stoch_eqs.noise_equations, stoch_eqs.operator_equations, stoch_eqs.states, stoch_eqs.operators, stoch_eqs.hamiltonian,stoch_eqs.jumps, stoch_eqs.jumps_dagger, stoch_eqs.rates, stoch_eqs.iv, stoch_eqs.varmap, 2)
test_mf_scaled = scale(test_mf)

for (lhs, rhs) in zip(stoch_eqs_scaled.noise_equations, test_mf_scaled.equations)
    @test isequal(simplify(lhs.lhs - rhs.lhs),0)
    @test isequal(simplify(lhs.rhs - rhs.rhs),0)
end

end