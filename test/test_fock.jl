using Qumulants
using SymbolicUtils
SymbolicUtils.show_simplified[]=false
using Test

# Test
hf = FockSpace(:c)

a = Destroy(hf,:a)
ad = a'

# Test conversion to Symbolics
a_sym = Qumulants._to_symbolic(a)
@test Qumulants.acts_on(a_sym)==1
@test Qumulants._to_qumulants(a_sym)==a

b = Destroy(hf,:b)
@test Qumulants._to_symbolic(a) != Qumulants._to_symbolic(b)

@test a==ad'
@test simplify_operators(a)==a
@test simplify_operators(a+a)==2*a
@test simplify_operators(a*a') == 1+a'*a
@test simplify_operators(a*a' + 1) == 2 + a'*a


# @syms x
# ex = x*a+x*a
# @test simplify_operators(ex)==simplify_operators((2*a)*x)
# @test isa(2*x*a, OperatorTerm) # TODO: Variadic methods -- use own Constants type

# Single mode
ωc = 0.1313
H = ωc*a'*a
da = simplify_operators(1.0im*(H*a - a*H))
@test da == -1.0im*ωc*a
