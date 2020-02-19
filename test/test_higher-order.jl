using Qumulants
using Test

@testset "higher-order" begin

# Test single mode
a = Destroy(:a)

H = a + a' + a'*a'*a*a # Kerr-nonlinearity
J = [a]

op = a'*a'*a*a
eq = heisenberg(op,H)
eq2 = heisenberg(op,H,J)

# Test single-atom laser
σ(i,j) = Transition(:σ,i,j,2)
a = Destroy(:a) ⊗ Identity()
s = Identity() ⊗ σ(1,2)
sz = simplify_operators(2*s'*s - one(s))
sps = simplify_operators(s'*s)
ωc, ωa, g, γ, κ, ν = (0.1, 3, 0.5, 0.25, 1.0, 4)
H = ωc*a'*a + ωa*s'*s + g*(a'*s + a*s')
J = [sqrt(κ)*a,sqrt(γ)*s,sqrt(ν)*s']

op = a'*a*a*s'

eq = heisenberg(op,H)
average(eq,4)

end # testset
