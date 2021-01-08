using Qumulants
using OrdinaryDiffEq
using Test

@testset "scaling" begin

M = 2

# Define parameters
@parameters Δ g γ κ ν

# Define hilbert space
hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atom, i),(:g,:e)) for i=1:M]
h = ⊗(hf, ha...)

# Define the fundamental operators
a = Destroy(h,:a)
σ(i,j,k) = Transition(h,Symbol(:σ_, k),i, j, k+1)

# Hamiltonian
H = Δ*a'*a + g*sum(a'*σ(:g,:e,i) + a*σ(:e,:g,i) for i=1:M)

# Collapse operators
J = [a;[σ(:g,:e,i) for i=1:M]; [σ(:e,:g,i) for i=1:M]]
rates = [κ;[γ for i=1:M];[ν for i=1:M]]

# Derive equation for average photon number
he_n = average(heisenberg(a'*a,H,J;rates=rates), M)

@parameters N
he_scale = scale(he_n, [2:M+1;], 1, N)
@test he_scale.rhs[1] == simplify_constants(-1.0κ*average(a'*a) + (-1.0im)*N*g*average(a'*σ(:g,:e,1)) + 1.0im*N*g*average(a*σ(:e,:g,1)))

# Custom filter function -- include only phase-invaraint terms
ϕ(x) = 0
ϕ(x::Destroy) = -1
ϕ(x::Create) = 1
function ϕ(t::Transition)
    if t.i < t.j
        1
    elseif t.i > t.j
        -1
    else
        0
    end
end
ϕ(avg::Average) = ϕ(avg.operator)
function ϕ(t::OperatorTerm)
    @assert t.f === (*)
    p = 0
    for arg in t.arguments
        p += ϕ(arg)
    end
    return p
end
phase_invariant(x) = iszero(ϕ(x))

# Complete equations
he = complete(he_scale;order=M,filter_func=phase_invariant);

find_missing(he.rhs, he.lhs)

# Correlation function
c = CorrelationFunction(a', a, he; steady_state=true, filter_func=phase_invariant)

# Numerical solution
ps = (Δ, g, γ, κ, ν, N)
f = generate_ode(he,ps)
u0 = zeros(ComplexF64, length(he))
p0 = (0, 1.5, 0.25, 1, 4, 50)
prob = ODEProblem(f,u0,(0.0,10.0),p0)
sol = solve(prob,RK4());

# Spectrum
S = Spectrum(c,ps)
ω = range(-pi, pi, length=501)
s = S(ω,sol.u[end],p0)

end # testset
