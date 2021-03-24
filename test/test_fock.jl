using Qumulants
using SymbolicUtils
using Test

@testset "fock" begin

# Test
hf = FockSpace(:c)

a = Destroy(hf,:a)
ad = a'

@test !isequal(hash(a), hash(ad))
b = Destroy(hf,:b)
@test !isequal(hash(a), hash(b))

@test isequal(a,ad')
@test isequal(simplify(a+a),2*a)
@test isequal(simplify(a/2 + 0.5*a),a)
@test isequal(a*a' , 1+a'*a)
@test isequal(simplify(a*a' + 1) , 2 + a'*a)

@test isequal(simplify(-1*(a'*a + 1)*a + a) , -1*a'*a^2)
@test isequal(simplify(a'*a*a - a*a'*a) , -1*a)

# Single mode
ωc = 0.1313
H = ωc*a'*a
da = simplify(1.0im*(H*a - a*H))
@test isequal(da , (0.0-1.0im)*ωc*a)

end # testset
