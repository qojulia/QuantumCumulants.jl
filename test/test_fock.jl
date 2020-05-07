using Qumulants
using SymbolicUtils
using Test

# Test
hf = FockSpace(:c)

a = Destroy(hf,:a)
ad = a'

@test a==ad'
@test simplify_operators(a)==a
@test simplify_operators(a+a)==2*a

@syms x
ex = x*a+x*a
@test simplify_operators(ex)==simplify_operators((2*a)*x)
# @test isa(2*x*a, OperatorTerm) # TODO: Variadic methods -- use own Constants type

tmp = Qumulants._to_symbolic(ex)

ex2 = -1*a*a'

ex3 = (a + a)*(a + a') + one(a)

# Single mode
ωc = 1.0
H = ωc*a'*a
da = simplify_operators(1.0im*(H*a - a*H))

tmp = Qumulants._to_symbolic(da)
