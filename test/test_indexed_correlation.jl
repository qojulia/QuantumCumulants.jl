using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics
using OrdinaryDiffEq
using SteadyStateDiffEq
using ModelingToolkit

const qc = QuantumCumulants

@testset "indexed_correlation" begin

order = 2 #order of the cumulant expansion
@cnumbers N Δ g κ Γ R ν M

# Hilbertspace
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom,2)

h = hc ⊗ ha

# Indices and Operators

k = Index(h,:k,N,ha)
l = Index(h,:l,N,ha)
m = Index(h,:m,N,ha)
n = Index(h,:n,N,ha)

extra_indices = [n]

@qnumbers a::Destroy(h)

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)

# Define the Hamiltonian
H = -Δ*a'a + g*(Σ(a'*σ(1,2,k),k) + Σ(a*σ(2,1,k),k))

J = [a,σ(1,2,l),σ(2,1,l),σ(2,2,l)]
rates = [κ, Γ, R, ν]
ops = [a'*a,σ(2,2,m)]
eqs = meanfield(ops,H,J;rates=rates,order=order)

@test length(eqs) == 2

φ(x::Average) = φ(x.arguments[1])
φ(::Destroy) = -1
φ(::Create) =1
φ(x::QTerm) = sum(map(φ, x.args_nc))
φ(x::Transition) = x.i - x.j
φ(x::IndexedOperator) = φ(x.op)
φ(x::SingleSum) = φ(x.term)
φ(x::AvgSums) = φ(arguments(x))
phase_invariant(x) = iszero(φ(x))

eqs_c = qc.complete(eqs;filter_func=phase_invariant);


eqs_sc1 = scale(eqs_c)

@test length(eqs_sc1) == 4
@test length(eqs_c) == 4


# define the ODE System and Problem
@named sys = ODESystem(eqs_sc1)

# Initial state
u0 = zeros(ComplexF64, length(eqs_sc1))
# System parameters
N_ = 2e5
Γ_ = 1.0 #Γ=1mHz
Δ_ = 2500Γ_ #Δ=2.5Hz
g_ = 1000Γ_ #g=1Hz
κ_ = 5e6*Γ_ #κ=5kHz
R_ = 1000Γ_ #R=1Hz
ν_ = 1000Γ_ #ν=1Hz

ps = [N, Δ, g, κ, Γ, R, ν]
p0 = [N_, Δ_, g_, κ_, Γ_, R_, ν_]

prob = ODEProblem(sys,u0,(0.0, 1.0/50Γ_), ps.=>p0);

# Solve the Problem
sol = solve(prob,Tsit5(),maxiters=1e7)

avrgSum = arguments(arguments(eqs_c[1].rhs)[1])[2]

split0 = split_sums(avrgSum,l,15)

@test isequal(split0,avrgSum)

split1 = split_sums(avrgSum,k,4)

@test split1 isa SymbolicUtils.BasicSymbolic && operation(split1) === *
@test !isequal(split1,avrgSum)
@test isequal(4,arguments(split1)[1])
@test isequal(N/4,arguments(split1)[2].metadata[qc.IndexedAverageSum].sum_index.range)

split2 = split_sums(avrgSum,k,M)

@test !isequal(split2,split1)
@test !isequal(split2,avrgSum)
@test isequal(M,arguments(split2)[1])
@test isequal(N/M,arguments(split2)[2].metadata[qc.IndexedAverageSum].sum_index.range)

# avrgSum2 = arguments(arguments(eqs_c[3].rhs)[1])[2] # order of arguments changed!
args_ = (arguments(eqs_c[3].rhs))
it_ = findfirst(x->typeof(arguments(x)[2]) == SymbolicUtils.BasicSymbolic{IndexedAverageSum}, args_)
avrgSum2 = arguments(arguments(eqs_c[3].rhs)[1it_])[2]

_split0 = split_sums(avrgSum2,l,15)

@test isequal(avrgSum2,_split0)
@test isequal(avrgSum2,split_sums(avrgSum2,k,1))

_split1 = split_sums(avrgSum2,k,3)

@test !isequal(_split1,avrgSum2)
@test _split1 isa SymbolicUtils.BasicSymbolic && operation(_split1) === +
@test isequal(2,length(arguments(_split1)))

op1 = a'
op2 = a

corr = IndexedCorrelationFunction(a', a, eqs_c; steady_state=true, filter_func=phase_invariant,extra_indices=extra_indices);
corr = scale(corr)

corr_nss = IndexedCorrelationFunction(a', a, eqs_c; steady_state=false, filter_func=phase_invariant);
corr_nss_sc = scale(corr_nss)

op1 = σ(1,2,m)
op2 = σ(2,1,n)

corr2 = CorrelationFunction(op1, op2, eqs_c; filter_func = phase_invariant) 
corr2_sc = scale(corr2);

p0_c2 = correlation_p0(corr2_sc, sol.u[end],ps.=>p0)
u0_c2 = correlation_u0(corr2_sc, sol.u[end])
@named csys2 = ODESystem(corr2_sc) # 5 equations, 8 parameters

# prob_c2 = ODEProblem(csys2,u0_c2,(0.0,5.0),p0_c2)
prob_c2 = ODEProblem(csys2,u0_c2,(0.0,0.05),p0_c2)
sol_c2 = solve(prob_c2,Tsit5();saveat=0.001,maxiters=1e8);

@test sol_c2.retcode == SciMLBase.ReturnCode.Success

# Test system with no filter func to evaluate corr function

eqs_c_2 = complete(eqs);
eqs_ev = evaluate(eqs_c_2;limits=(N=>3))

u0_ev = zeros(ComplexF64, length(eqs_ev))

@named sys_ev = ODESystem(eqs_ev)
# prob_ev = ODEProblem(sys_ev,u0_ev,(0.0, 1.0/50Γ_), ps.=>p0);
prob_ev = ODEProblem(sys_ev,u0_ev,(0.0, 0.1/50Γ_), ps.=>p0);
sol_ev = solve(prob_ev,Tsit5(),maxiters=1e7)
@test sol_ev.retcode == SciMLBase.ReturnCode.Success

corr3 = CorrelationFunction(op1, op2, eqs_c_2)
corr3_ev = evaluate(corr3,1,2; limits=(N=>3));

p0_c3 = correlation_p0(corr3_ev, sol_ev.u[end],ps.=>p0)
u0_c3 = correlation_u0(corr3_ev, sol_ev.u[end])
@named csys3 = ODESystem(corr3_ev) 

# prob_c3 = ODEProblem(csys3,u0_c3,(0.0,5.0),p0_c3)
prob_c3 = ODEProblem(csys3,u0_c3,(0.0,0.05),p0_c3)
sol_c3 = solve(prob_c3,Tsit5();saveat=0.001,maxiters=1e8);

@test sol_c3.retcode == SciMLBase.ReturnCode.Success

ps = [N, Δ, g, κ, Γ, R, ν]

@test length(corr.de0) == 4
@test length(corr.de) == 2

S = Spectrum(corr, ps);

# prob_ss = SteadyStateProblem(prob)
# sol_ss = solve(prob_ss, DynamicSS(Tsit5()),
#     reltol=1e-14, abstol=1e-14, maxiters=5e7);
sol_ss = solve(prob, Tsit5(), save_everystep=false, save_on=false, save_start=false)

p0_c = correlation_p0(corr_nss_sc, sol.u[end],ps.=>p0)
u0_c = correlation_u0(corr_nss_sc, sol.u[end])

@named csys = ODESystem(corr_nss_sc) 

# prob_c4 = ODEProblem(csys,u0_c,(0.0,10.0),p0_c);
prob_c4 = ODEProblem(csys,u0_c,(0.0,0.05),p0_c);
sol_c = solve(prob_c4,Tsit5();saveat=0.001,maxiters=1e8); 

@test sol_c.retcode == SciMLBase.ReturnCode.Success

limits = Dict{SymbolicUtils.BasicSymbolic,Int64}(N=>5)
evals = evaluate(eqs_c;limits=limits)

corr4 = CorrelationFunction(op1, op2, eqs_c; filter_func = phase_invariant, steady_state=true) 
corr4_sc = scale(corr4);

S4_sc = Spectrum(corr4_sc, ps)
ω = [-10:0.01:10;]Γ_
S4_sc(ω,sol_ss.u[end],p0)
S4_sc(ω,sol.u[end],p0)

@test length(evals) == 21
@test length(unique(evals.states)) == length(evals.states)

@cnumbers N_

i = Index(h,:i,N_,ha)
j = Index(h,:j,N_,ha)

Γ_ij = DoubleIndexedVariable(:Γ,i,j)
Ω_ij = DoubleIndexedVariable(:Ω,i,j;identical=false)
gi = IndexedVariable(:g,i)

ΓMatrix = [1.0 1.0
            1.0 1.0]

ΩMatrix = [0.0 0.5
            0.5 0.0]

gn = [1.0,2.0]
κn = 3.0

ps = [Γ_ij,Ω_ij,gi,κ]
p0 = [ΓMatrix,ΩMatrix,gn,κn]

limits = Dict{SymbolicUtils.BasicSymbolic,Int64}(N_=>2)
valmap = value_map(ps,p0;limits=limits)

map_ = Dict{SymbolicUtils.BasicSymbolic,ComplexF64}()
push!(map_,qc.DoubleNumberedVariable(:Γ,1,1)=>1.0)
push!(map_,qc.DoubleNumberedVariable(:Γ,2,1)=>1.0)
push!(map_,qc.DoubleNumberedVariable(:Γ,1,2)=>1.0)
push!(map_,qc.DoubleNumberedVariable(:Γ,2,2)=>1.0)
push!(map_,qc.DoubleNumberedVariable(:Ω,2,2)=>0.0)
push!(map_,qc.DoubleNumberedVariable(:Ω,1,2)=>0.5)
push!(map_,qc.DoubleNumberedVariable(:Ω,1,1)=>0.0)
push!(map_,qc.DoubleNumberedVariable(:Ω,2,1)=>0.5)
push!(map_,qc.SingleNumberedVariable(:g,1)=>1.0)
push!(map_,qc.SingleNumberedVariable(:g,2)=>2.0)
push!(map_,κ=>3.0)



@test length(valmap) == 11
@test length(valmap) == length(map_)
for arg in valmap
    @test arg in map_
end

@test isequal(qc.scale_term(σ(1,2,k)), σn(1,2,1))
@test isequal(qc.scale_term(σ(2,1,k)*σ(1,2,j)), σn(2,1,1)*σn(1,2,2))

sum_ = average(Σ(σ(2,1,i)*σ(1,2,j),i))

split = split_sums(sum_,i,3)

@test split isa SymbolicUtils.BasicSymbolic && operation(split) === +
@test length(arguments(split)) == 3

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom,2)

h = hc ⊗ ha

k = Index(h,:k,N,ha)
l = Index(h,:l,N,ha)

m = Index(h,:m,N,hc)
n = Index(h,:n,N,hc)

order = 2

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
ai(k) = IndexedOperator(Destroy(h,:a),k)

H_2 = -Δ*∑(ai(m)'ai(m),m) + g*(∑(Σ(ai(m)'*σ(1,2,k),k),m) + ∑(Σ(ai(m)*σ(2,1,k),k),m))

J_2 = [ai(m),σ(1,2,k),σ(2,1,k),σ(2,2,k)]
rates_2 = [κ, Γ, R, ν]
ops_2 = [ai(n)'*ai(n),σ(2,2,l)]

q = Index(h,:q,N,ha)
r = Index(h,:r,N,hc)

extra_indices = [q,r]
eqs_2 = meanfield(ops_2,H_2,J_2;rates=rates_2,order=order)

eqs_s1 = scale(eqs_2;h=2)
inds_s1 = qc.get_indices_equations(eqs_s1)

for ind in inds_s1
    @test isequal(ind.aon,1)
end

eqs_s2 = scale(eqs_2;h=1)
inds_s2 = qc.get_indices_equations(eqs_s2)

for ind in inds_s2
    @test isequal(ind.aon,2)
end

sum1 = ∑(σ(2,1,i),i)
op1 = a'
op2 = a
h1 = qc.hilbert(op1)
h2 = qc._new_hilbert(qc.hilbert(op2), acts_on(op2))
h_ = h1⊗h2
new_sum = qc._new_operator(sum1,h_)
@test isequal(qc.hilbert(new_sum),h_)

# Test that checks for correct spectrum, when using multiple indexed Hilbert spaces
# Hamiltonian is based on the indexed superrad laser
    hc = FockSpace(:cavity)
    hA = NLevelSpace(:atomA,2)
    hB = NLevelSpace(:atomB,2)
    h = hc ⊗ hA ⊗ hB
    # operators
    @qnumbers a::Destroy(h)
    σA(α,β,i) = IndexedOperator(Transition(h,Symbol("σ_a"), α, β, 2),i)
    σB(α,β,i) = IndexedOperator(Transition(h,Symbol("σ_b"), α, β, 3),i)
    @cnumbers N Δ κ Γ R ν
    g(i) = IndexedVariable(:g, i)
    i = Index(h,:i,N,hA)
    j = Index(h,:j,N,hA)
    x = Index(h,:x,N,hB)
    y = Index(h,:y,N,hB)
    # Hamiltonian
    H = -Δ*a'a + Σ(g(i)*( a'*σA(1,2,i) + a*σA(2,1,i) ),i) + Σ(σA(1,2,i)*σB(2,1,x) + σA(2,1,i)*σB(1,2,x),i,x) + Σ(g(x)*( a'*σB(1,2,x) + a*σB(2,1,x) ),x)

    # Jump operators wth corresponding rates
    J = [a, σA(1,2,i), σA(2,1,i), σA(2,2,i), σB(1,2,x), σB(2,1,x), σB(2,2,x)]
    rates = [κ, Γ, R, ν, Γ, R, ν]
    # Derive equations
    ops = [a'*a, σA(2,2,j)]
    eqs = meanfield(ops,H,J;rates=rates,order=2);
    eqs_c = complete(eqs);
    eqs_sc = scale(eqs_c);
    @named sys = ODESystem(eqs_sc)
    # Initial state
    u0 = zeros(ComplexF64, length(eqs_sc))
    # System parameters
    N_ = 2e5
    Γ_ = 1.0 #Γ=1mHz
    Δ_ = 2500Γ_ #Δ=2.5Hz
    g_ = 1000Γ_ #g=1Hz
    κ_ = 5e6*Γ_ #κ=5kHz
    R_ = 1000Γ_ #R=1Hz
    ν_ = 1000Γ_ #ν=1Hz
    ps = [N, Δ, g(1), κ, Γ, R, ν]
    p0 = [N_, Δ_, g_, κ_, Γ_, R_, ν_]
    # prob = ODEProblem(sys,u0,(0.0, 1.0/50Γ_), ps.=>p0)
    prob = ODEProblem(sys,u0,(0.0, 0.1/50Γ_), ps.=>p0)
    # Solve the numeric problem
    sol = solve(prob,Tsit5(),maxiters=1e7);
    b = Index(h,:b,N,hA)
    lims = (N=>2)
    op1 = σA(1,2,j)
    op2 = σA(2,1,b)
    corr = CorrelationFunction(op1, op2, eqs_c; steady_state=true) 
    corr_sc = scale(corr;limits=lims);
    S_sc = Spectrum(corr_sc, ps)
    ω = [-10:0.01:10;]Γ_
    S_sc(ω,sol.u[end],p0);
    # prob_ss = SteadyStateProblem(prob)
    # sol_ss = solve(prob_ss, DynamicSS(Tsit5()),
    #     reltol=1e-14, abstol=1e-14, maxiters=5e7);

    @test sol.retcode == SciMLBase.ReturnCode.Success
    # @test sol_ss.retcode == SciMLBase.ReturnCode.Success
# end of testcase

end
