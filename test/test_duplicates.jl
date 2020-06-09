using Qumulants
using OrdinaryDiffEq
using Test
using Random

@testset "duplicates" begin

Random.seed!(0)

# Hilbert space
N = 4
hf = FockSpace(:field)
ha = NLevelSpace(:atom,(:g,:e))
h_tot = hf⊗ha

# Operators
a = Destroy(h_tot,:a)
σ(i,j,k) = k isa Index ? Transition(h_tot, :σ, i, j; index=k) : Transition(h_tot, :σ, i, j; index=I[k])
σ(k) = σ(:g,:e,k)

# Indices
I = IndexSet(:I, 1, N)

# Test simple system
H = +((σ(i)'*σ(i) for i=I)...) + +((σ(i) + σ(i)' for i=I)...)
ops = [σ(I[1]), σ(:e,:e,I[1])]
he = heisenberg(ops, H)
he_dupl = build_duplicates(he)

# Compare to direct derivation
ops_comp = [[σ(i) for i=I]; [σ(:e,:e,i) for i=I];]
he_comp  = heisenberg(ops_comp, H)

@test length(he_dupl)==length(he_comp)
for i=1:length(he_dupl)
    j = findfirst(isequal(he_dupl.lhs[i]), he_comp.lhs)
    @test he_dupl.rhs[i] == he_comp.rhs[j]
end

# Test with indexed parameters
@parameters Δ Ω
H = +((Δ[i]*σ(i)'*σ(i) for i=I)...) + +((Ω[i]*(σ(i) + σ(i)') for i=I)...)
he = heisenberg(ops, H)
he_dupl = build_duplicates(he)
he_comp = heisenberg(ops_comp, H)

@test length(he_dupl)==length(he_comp)
for i=1:length(he_dupl)
    j = findfirst(isequal(he_dupl.lhs[i]), he_comp.lhs)
    @test he_dupl.rhs[i] == he_comp.rhs[j]
end

# Many-atom laser equations
@parameters Δ g κ γ ν
H = +((Δ[i]*σ(i)'*σ(i) for i=I)...) + +((g[i]*(a'*σ(i) + a*σ(i)') for i=I)...)
J = [a;[σ(i) for i=I];[σ(i)' for i=I]]
rates = [κ;[γ[i] for i=I];[ν[i] for i=I]]

ops = [a'*a,a'*σ(I[1]),σ(:e,:e,I[1]),σ(:e,:g,1)*σ(:g,:e,2)]
he = heisenberg(ops,H,J;rates=rates)
he_dupl = build_duplicates(he)

ops_comp = [a'*a; [a'*σ(i) for i=I]; [σ(:e,:e,i) for i=I]]
for j=1:N
    for k=j+1:N
        if j==1 && k==2
            push!(ops_comp, σ(:e,:g,j)*σ(:g,:e,k))
        else
            push!(ops_comp, σ(:g,:e,j)*σ(:e,:g,k))
        end
    end
end

@test length(ops_comp)==length(he_dupl)
@test all(op in he_dupl.lhs for op in ops_comp)
@test all(op in ops_comp for op in he_dupl.lhs)

he_comp = heisenberg(ops_comp,H,J;rates=rates)
for i=1:length(he_dupl)
    j = findfirst(isequal(he_dupl.lhs[i]), he_comp.lhs)
    @test simplify_operators(he_dupl.rhs[i]) == he_comp.rhs[j]
end

# Test averaging
he_avg = average(he, 2)
he_dupl_avg = build_duplicates(he_avg)
he_comp_avg = average(he_comp, 2)
@test length(he_dupl_avg)==length(he_comp_avg)
for i=1:length(he_dupl_avg)
    j = findfirst(isequal(he_dupl_avg.lhs[i]), he_comp_avg.lhs)
    # TODO bug in sorting?
    # @test simplify_constants(he_dupl_avg.rhs[i]) == simplify_constants(he_comp_avg.rhs[j])
end

# Test numerical solution for many-atom laser
ps_dupl = [κ;[Δ[i] for i=I];[g[i] for i=I];[γ[i] for i=I];[ν[i] for i=I]]
missed = find_missing(he_dupl_avg;ps=ps_dupl)
@test !any(m isa Parameter for m in missed)
subs = Dict(missed .=> 0) # phase invariance
he_nophase = simplify_constants(substitute(he_dupl_avg, subs))
f_dupl = generate_ode(he_nophase, ps_dupl)

u0 = zeros(ComplexF64, length(he_nophase))
p_numeric = abs.(ones(length(ps_dupl)) .+ randn(length(ps_dupl)))
prob = ODEProblem(f_dupl, u0, (0.0, 10.0), p_numeric)
sol_dupl = solve(prob, RK4())

n_dupl = getindex.(sol_dupl.u, findfirst(isequal(a'*a), he_dupl.lhs))
pops_dupl = [getindex.(sol_dupl.u, findfirst(isequal(σ(:e,:e,i)), he_dupl.lhs)) for i=1:N]
@test all(abs.(imag.(n_dupl)) .< 1e-16)
for i=1:N
    @test all(abs.(imag.(pops_dupl[i])) .< 1e-16)
    @test all(1.0 .>= real.(pops_dupl[i]) .>= 0.0)
end


# "Brute-force" implementation
@parameters κ
Δ = [parameters([Symbol(:Δ, i) for i=1:N]...)...]
g = [parameters([Symbol(:g, i) for i=1:N]...)...]
γ = [parameters([Symbol(:γ, i) for i=1:N]...)...]
ν = [parameters([Symbol(:ν, i) for i=1:N]...)...]
hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
h = hf ⊗ ha

a = Destroy(h,:a,1)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k+1)
H = sum(Δ[i]*σ(2,2,i) for i=1:N) + sum(g[i]*(a'*σ(1,2,i) + a*σ(2,1,i)) for i=1:N)
J = [a;[σ(1,2,i) for i=1:N];[σ(2,1,i) for i=1:N]]
rates = [κ;γ;ν]
ops = [a'*a; [a'*σ(1,2,i) for i=1:N]; [σ(2,2,i) for i=1:N]]
for j=1:N
    for k=j+1:N
        push!(ops, σ(2,1,j)*σ(1,2,k))
    end
end
he = heisenberg(ops, H, J; rates=rates)
he_avg = average(he,2)
ps_brute = [κ;Δ;g;γ;ν]
missed = find_missing(he_avg;ps=ps_brute)
@test !any(m isa Parameter for m in missed)
subs = Dict(missed .=> 0)
he_nophase = simplify_constants(substitute(he_avg, subs))
f_brute = generate_ode(he_nophase, ps_brute)

prob = ODEProblem(f_brute, u0, (0.0, 10.0), p_numeric)
sol_brute = solve(prob, RK4())

@test length(sol_dupl)==length(sol_brute)
n_brute = getindex.(sol_brute.u, findfirst(isequal(a'*a), he.lhs))
pops_brute = [getindex.(sol_brute.u, findfirst(isequal(σ(2,2,i)), he.lhs)) for i=1:N]
@test all(abs.(imag.(n_brute)) .< 1e-16)
for i=1:N
    @test all(abs.(imag.(pops_brute[i])) .< 1e-16)
    @test all(1.0 .>= real.(pops_brute[i]) .>= 0.0)
end

@test isapprox(n_brute,n_dupl)
@test isapprox(pops_brute,pops_dupl)

end # testset
