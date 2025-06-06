using QuantumCumulants
using SymbolicUtils
using Test

@testset "fock" begin

# Test
hf = FockSpace(:c)

a = Destroy(hf,:a)
P0 = GroundStateProjection(hf, :a)
ad = a'

@test !isequal(hash(a), hash(ad))
b = Destroy(hf,:b)
@test !isequal(hash(a), hash(b))
@test !isequal(hash(a), hash(P0))

@test isequal(a,ad')
@test isequal(simplify(a+a),2*a)
@test isequal(simplify(a/2 + 0.5*a),a)
@test isequal(a*a' , 1+a'*a)
@test isequal(simplify(a*a' + 1) , 2 + a'*a)
@test isequal(simplify(P0 + P0) , 2*P0)
@test isequal(simplify(P0^2) , P0)
@test isequal(simplify(a'*a*(a')^3*P0^2) , 3*(a')^3*P0)

@test isequal(simplify(-1*(a'*a + 1)*a + a) , -1*a'*a^2)
@test isequal(simplify(a'*a*a - a*a'*a) , -1*a)
@test isequal(simplify(P0 - P0) , 0)

# Single mode
ωc = 0.1313
H = ωc*a'*a
da = simplify(1.0im*(H*a - a*H))
@test isequal(da , (0.0-1.0im)*ωc*a)

# Test substitute by numbers
@syms x
@test iszero(substitute(x*a, Dict(x=>0)))
@test isequal(substitute(x*a, Dict(x=>1)), a)
@test iszero(substitute(x*a, Dict(a=>0)))
@test iszero(substitute(x*P0, Dict(x=>0)))
@test isequal(substitute(x*P0, Dict(x=>1)), P0)
@test iszero(substitute(x*P0, Dict(P0=>0)))

# Test substitute by syms
@syms y
@test isequal(substitute(x*a, Dict(x=>y)), y*a)
@test isequal(substitute(x*a, Dict(a=>y)), x*y)
@test isequal(substitute(x*(a+a'), Dict(x => y)), y*(a + a'))
@test isequal(substitute(x*P0, Dict(x=>y)), y*P0)
@test isequal(substitute(x*P0, Dict(P0=>y)), x*y)
@test isequal(substitute(x*(P0+a'*P0), Dict(x => y)), y*(P0 + a'*P0))

end # testset
