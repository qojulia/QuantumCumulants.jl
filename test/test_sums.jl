using Qumulants
using Test
using SymPy
using MacroTools

δ(i,j) = SymPy.sympy.functions.special.tensor_functions.KroneckerDelta(sympify(i),sympify(j))
#
# # @testset "sums" begin

# # Test single mode
# a = Destroy(:a)
# σ(i,j) = Transition(:σ,i,j,2)
# i = Index(:i,1,:N)
# j = Index(:j,1,:N)
#
# S = Sum(a[i],i)
# @test S==Sum(a[i])
# @test S==simplify_operators(S)
#
# S = Sum(δ(i,j)*a[i],i)
# @test simplify_operators(S)==a[j]
# S = Sum(δ(i,j)*a[j],i)
# @test simplify_operators(S)==a[j]
# S = Sum(δ(i,j)*a[j],j)
# @test simplify_operators(S)==a[i]
#
# g = symbols("g",real=true)
# S = Sum(δ(i,j)*g[i]*a[i]*a[j],i)
# @test simplify_operators(S)==g[j]*a[j]^2
#
# S = Sum(δ(i,j)*(1-δ(i,j))*g[i]*a[i]*a[j],i)
# @test iszero(simplify_operators(S))
# S = Sum((g[j] + 3.333 + δ(i,j))*g[i]*a[j],j)
# @test simplify_operators(S)==g[i]*a[i] + Sum((g[j]*g[i]+3.333*g[i])*a[j],j)
# S = Sum((g[j] + 2 + δ(i,j))*g[i]*a[j],i)
# @test simplify_operators(S)== g[j]*a[j] + Sum((g[i]*g[j]+2*g[i])*a[j],i)
#
# S = Sum(g[j]*δ(i,j)*a'⊗σ(2,2)[i], j)
# @test simplify_operators(S)==g[i]*a'⊗σ(2,2)[i]
#
# # Implement the Tavis-Cummings Hamiltonian
# summand = g[i]*(a'⊗σ(1,2)[i] + a⊗σ(2,1)[i])
#
# @test Qumulants.swap_index(g[i]*a[i],i,j) == g[j]*a[j]
# @test Qumulants.swap_index(-g[i]*a[i],i,j) == -g[j]*a[j]
#
# H_TC = Sum(summand,i)
# @test H_TC == Sum(g[i]*a'⊗σ(1,2)[i],i) + Sum(g[i]*a⊗σ(2,1)[i],i)
# @test iszero(simplify_operators(H_TC - H_TC))
#
# S1 = Sum(g[i]*a[i], i)
# S2 = Sum(g[j]*a[j], j)
# @test iszero(simplify_operators(S1 - S2))
#
# # Test heisenberg
# ops = [a⊗Identity(),Identity()⊗σ(1,2)[j]]
#
# tmp = g[j]*a⊗σ(2,1)[j]
#
# S = Sum(tmp,j)
# tmp2 = 1.0im*S*ops[2]
# tmp3 = ops[2]*S
#
# da = heisenberg(ops[1],H_TC)
# ds = heisenberg(ops[2],H_TC)
#
# ops = [a⊗Identity(),Identity()⊗σ(1,2)[j]]
# he = heisenberg(ops,H_TC)
#
# tmp = (1 - δ(i,j))*g[i]*a ⊗(σ(2,1)[i])
# _, tmp2 = Qumulants.resolve_kdelta(tmp,i)
#
# # tmp = g[j]*a⊗(δ(sympify(i),sympify(j))*σ(1,2)[j])
#
# tmp = (ops[2]*S.args[1]).args[1]
# ci, ti = Qumulants.find_index(tmp)
#
# # end # testset
#
#
# ##################################################################
# #christoph tests
# J_test = [a⊗Identity(),Identity()⊗σ(1,2)[i]]#[a,σ(2,1)[j]]
# Jd_test = adjoint.(J_test)
# master_test3 = heisenberg([J_test[1], Identity()⊗σ(1,2)[j]], H_TC, J_test;Jdagger=Jd_test)
#
# κ = symbols("κ",real=true,positive=true)
# γ = symbols("γ",real=true,positive=true)
# ops = [a⊗Identity(),Identity()⊗σ(1,2)[j],a'*a⊗Identity(),a'⊗σ(1,2)[j]]
# me = heisenberg(ops[end], H_TC, J_test;Jdagger=Jd_test,rates=[κ,γ[i]])
# he = heisenberg(ops[end],H_TC)
#
# tmp = simplify_operators(1.0im*(Sum(g[i]*a⊗σ(2,1)[i],i)*(a'⊗σ(1,2)[j]) - (a'⊗σ(1,2)[j])*Sum(g[i]*a⊗σ(2,1)[i],i)))
# tmp2 = Qumulants.resolve_kdelta(simplify_operators(tmp.args[1].args[1]),i)[2]
#
#
# # Test second order indexed
# k = Index(:k,1,:N)
# he = heisenberg(Identity()⊗(σ(1,2)[j]*σ(2,1)[k]),H_TC)
#
# tmp = Sum((-δ(j,k)*g[i] + g[i])*a⊗(σ(2,1)[i]*σ(1,2)[j]*σ(2,1)[k]),i)
#
#
# # Test double sums
# σ(i,j) = Transition(:σ,i,j,2)
# i = Index(:i,1,:N)
# j = Index(:j,1,:N)
#
# Ω = symbols("Ω",real=true)
# Δ = symbols("Δ",real=true)
# γ = symbols("γ",real=true,positive=true)
#
# H0 = Sum(Δ[i]*σ(2,2)[i],i)
# H = Sum(Sum(Ω[i,j]*σ(2,1)[i]*σ(1,2)[j],i),j)
#
# r = Index(:r,1,:N)
# s = Index(:s,1,:N; neq=[r])
# ops = [σ(1,2)[r], σ(2,2)[r], σ(2,1)[r]*σ(1,2)[s]]
# he = heisenberg(ops,H)
#
# tmp = Sum(Sum(δ(i,r)*σ(2,1)[i]*σ(1,2)[j],i),j)
# tmp2 = Sum(simplify_operators(tmp.args[1]), j)
#
#
# J = [σ(1,2)[i]]
# rates = [γ[i,j]]
# # ops = [σ(2,1)[r]*σ(1,2)[s]]
# he = heisenberg(ops[1],H)
#
# tmp = he.rhs.args[2]
#
# he1 = heisenberg(σ(2,2)[r],simplify_operators(H))
#
# he2 = heisenberg(σ(1,2)[r],H)
# tmp1, tmp2 = he2.rhs.args[3], he2.rhs.args[5]
#
# me = heisenberg(ops[1], 0H, J; rates=rates)
#
# op1 = σ(2,1)[i]
# op2 = σ(1,2)[j]
# tmp = γ[i,j]*(Qumulants.commutator(σ(1,2)[r],op1)*op2)
# tmp2 = Sum(Sum((tmp),i), j).args[3]


# # heisenberg equation up to 2nd order
# a = Destroy(:a)
# σ(i,j) = Transition(:σ,i,j,2)
# i = Index(:i,1,:N)
# j = Index(:j,1,:N)
# k = Index(:k,1,:N;neq=[j])
#
# g = symbols("g",real=true)
# Δ = symbols("Δ",real=true)
# κ = symbols("κ",real=true,positive=true)
# γ = symbols("γ",real=true,positive=true)
#
# # Implement the Tavis-Cummings Hamiltonian
# H_int_i = g[i]*(a'⊗σ(1,2)[i] + a⊗σ(2,1)[i])
# H_TC = Δ*(a'*a)⊗Identity() + Sum(H_int_i,i)
# J = [a⊗Identity(),Identity()⊗σ(1,2)[i]];
# ops = [a⊗Identity(),Identity()⊗σ(1,2)[j],(a'*a)⊗Identity(),a'⊗σ(1,2)[j],
#     Identity()⊗σ(2,2)[j],Identity()⊗(σ(2,1)[j]*σ(1,2)[k])]
# me = heisenberg(ops,H_TC,J;rates=[κ,γ[i]])
# tmp = me.rhs[end].args[1]

# Define Parameters
g, Δ, κ, γ, ν = symbols("g Δ κ γ ν", real=true)

# Define Operators
a = Destroy(:a) ⊗ Identity()
σ(i,j,k) = Identity() ⊗ Transition(:σ,i,j,(:g,:e))[k]

# Define symbolic Indices
i = Index(:i,1,:n)
j = Index(:j,1,:n)
k = Index(:k,1,:n;neq=[j])

# Implement the Tavis-Cummings Hamiltonian and decay
H_TC = Δ*(a'*a) + Sum(g[i]*(a'*σ(:g,:e,i) + a*σ(:e,:g,i)),i)
J = [a,σ(:g,:e,i),σ(:e,:g,i)]

ops = [(a'*a),a'*σ(:g,:e,j),σ(:e,:e,j),σ(:e,:g,k)*σ(:g,:e,j)]
he = heisenberg(ops,H_TC,J;rates=[κ,γ[i],ν[i]])

he_avg = average(he,2)

# tmp = Qumulants.order_indexed(he_avg.rhs,ops,2)

n = symbols("n", integer=true)
p = (Δ,κ,g,γ,ν,n-1)
meta_f = build_ode(he_avg,p;set_unknowns_zero=true)


tmp3 = average(ops[end])
tmp2 = tmp3.__pyobject__.base[sympify(Qumulants.IndexOrder[6]),sympify(Qumulants.IndexOrder[4])]

tmp4 = meta_f[2].__pyobject__.replace(tmp2,SymPy.symbols("HERE"))

using Combinatorics
tmp2 = combinations(Qumulants.IndexOrder,2)
