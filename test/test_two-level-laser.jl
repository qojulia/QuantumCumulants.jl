using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using Test

@testset "two-level-laser" begin

N = 10

Δ = cnumbers((Symbol(:Δ_, i) for i=1:N)...)
g = cnumbers((Symbol(:g_, i) for i=1:N)...)
γ = cnumbers((Symbol(:γ_, i) for i=1:N)...)
ν = cnumbers((Symbol(:ν_, i) for i=1:N)...)
@cnumbers κ

h_cavity = FockSpace(:cavity)
h_atoms = [NLevelSpace(Symbol(:atom, i),(:g,:e)) for i=1:N]
h = tensor(h_cavity, h_atoms...)

a = Destroy(h, :a)
σ(i,j,k) = Transition(h,Symbol(:σ_,k),i,j,k+1)

ops = begin
    ops_ = [a'*a;[σ(:e,:e,k) for k=1:N];[a'*σ(:g,:e,k) for k=1:N]]
    for i=1:N
        for j=i+1:N
            push!(ops_, σ(:e,:g,i)*σ(:g,:e,j))
        end
    end
    ops_
end

n_eqs = div(N*(N-1),2) + 2N + 1

H = sum(Δ[i]*σ(:e,:e,i) for i=1:N) + sum(g[i]*(a'*σ(:g,:e,i) + a*σ(:e,:g,i)) for i=1:N)

J = [a;[σ(:g,:e,k) for k=1:N];[σ(:e,:g,k) for k=1:N]]
Jdagger = adjoint.(J)
rates = [κ,γ...,ν...]
he = heisenberg(ops, H, J; Jdagger=Jdagger, rates=rates, simplify=true, expand=true, order=2)

missed = find_missing(he)

subs = Dict(missed .=> 0)

he_nophase = substitute(he, subs)

@test isempty(find_missing(he_nophase))

eqs_mtk = equations(he_nophase)

sys = ODESystem(he_nophase)

u0 = zeros(ComplexF64, length(ops))
p0 = [κ => 1,
    (γ .=> 0.25 .* ones(N))...,
    (ν .=> 4 .* ones(N))...,
    (g .=> 1.5 .* ones(N))...,
    (Δ .=> ones(N))...,]

prob = ODEProblem(sys,u0,(0.0,10.0),p0,jac=true,sparse=true)

sol = solve(prob,RK4())

n = real.(sol[average(a'*a)])
pe = sol[average(σ(:e,:e,1))]


# Test with complete and custom filter
ϕ(::Destroy) = -1
ϕ(::Create) = 1
function ϕ(t::Transition)
    if t.i != t.j
        t.i == :e && return 1
        return -1
    else
        return 0
    end
end
function ϕ(q::QuantumCumulants.QMul)
    p = 0
    for arg ∈ q.args_nc
        p += ϕ(arg)
    end
    return p
end
ϕ(avg::Average) = ϕ(avg.arguments[1])
phase_invariant(x) = iszero(ϕ(x))

he_n = heisenberg(a'*a, H, J; rates=rates)
complete!(he_n;filter_func=phase_invariant)

@test length(he_n.equations) == length(ops)
@test isempty(find_missing(he_n))


sys_comp = ODESystem(he_n)
prob_comp = ODEProblem(sys_comp,u0,(0.0,10.0),p0)

sol_comp = solve(prob_comp,RK4())

@test getindex.(sol.u, 1) ≈ getindex.(sol_comp.u, 1)

end # testset
