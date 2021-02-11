using Qumulants
using OrdinaryDiffEq

M = 2
hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atom, i), (:g,:e)) for i=1:M]
h = tensor(hf, ha...)

@qnumbers a::Destroy(h)
σ(i,j,k) = Transition(h,Symbol(:σ_,k),i,j,k+1)

@cnumbers Δ g κ γ ν

H = Δ*a'*a + g*sum(a'*σ(:g,:e,i) + a*σ(:e,:g,i) for i=1:M)
J = [a;[σ(:g,:e,i) for i=1:M];[σ(:e,:g,i) for i=1:M]]
rates = [κ; [γ for i=1:M]; [ν for i=1:M]]

ops = [a'*a, a'*σ(:g,:e,1), σ(:e,:e,1), σ(:e,:g,1)*σ(:g,:e,2)]

he = heisenberg(ops, H, J; rates=rates)

# ϕ(x) = 0
# ϕ(::Destroy) = -1
# ϕ(::Create) = 1
# function ϕ(t::Transition)
#     if (t.i==:e && t.j==:g)
#         1
#     elseif (t.i==:g && t.j==:e)
#         -1
#     else
#         0
#     end
# end
# ϕ(avg::Average) = ϕ(avg.arguments[1])
# function ϕ(t::QTerm)
#     @assert t.f === (*)
#     p = 0
#     for arg in t.arguments
#         p += ϕ(arg)
#     end
#     return p
# end
# phase_invariant(x) = iszero(ϕ(x))
# missed = find_missing(he_avg)
# filter!(!phase_invariant, missed)
#
# he_nophase = substitute(he_avg, Dict(missed .=> 0))


# Scale
@cnumbers N
he_scaled = scale(he, [2,3], N)

ps = (Δ, g, γ, κ, ν, N)
f = generate_ode(he_scaled, ps)
p0 = (0, 1.5, 0.25, 1, 4, 7)
u0 = zeros(ComplexF64, length(he_scaled))
prob = ODEProblem(f, u0, (0.0, 50.0), p0)
sol = solve(prob, RK4())

ex = a*σ(:e,:g,1)*σ(:e,:e,2)
ex2 = a'*σ(:e,:e,1)*σ(:g,:e,2)
Qumulants.is_redundant(ex,[2,3])
ex_ = Qumulants._replace_redundant(ex,[3,2],[:a,:σ_1,:σ_2])
