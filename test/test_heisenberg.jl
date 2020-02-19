using Qumulants
using Test

@testset "heisenberg" begin

# Test single mode
a = Destroy(:a)


tmp = a*a'
tmp2 = Qumulants.apply_comms(tmp)
@test tmp2 == a'*a + one(a)

# Without SymPy
ω, κ = (10.0, 0.1)
H = ω*a'*a
J = [sqrt(κ)*a]

tmp = simplify_operators(a*H)
tmp2 = Qumulants.apply_comms(tmp)
@test tmp2 == H*a + ω*a

da = 1.0im*(H*a - a*H)
da_sim = Qumulants.apply_comms(da)
@test heisenberg(a,H).rhs == da_sim == -1.0im*ω*a
da_qle = simplify_operators(1.0im*(H*a - a*H) + sum(j'*a*j - 0.5*(j'*j*a + a*j'*j) for j=J))
da_qle_sim = Qumulants.apply_comms(da_qle)
@test heisenberg(a,H,J).rhs == da_qle_sim == -(1.0im*ω + 0.5κ)*a

# Test composite system
b = Destroy(:b)
id = Identity()

# With zero coupling -- should result in the same as before (⊗𝟙)
ωa, ωb, κa, κb, g = (10.0,6.0,1.0,1.1,0.0)

H = ωa*(a'*a)⊗id + ωb*id⊗(b'*b) + g*(a'⊗b + a⊗b')
J = [sqrt(κa)*a⊗id,sqrt(κb)*id⊗b]

a_ = a⊗id
da = simplify_operators(1.0im*(H*a_ - a_*H))
da_sim = Qumulants.apply_comms(da)
da_qle = simplify_operators(1.0im*(H*a_ - a_*H) + sum(j'*a_*j - 0.5*(j'*j*a_ + a_*j'*j) for j=J))
da_qle_sim = Qumulants.apply_comms(da_qle)
@test da_qle_sim == heisenberg(a_,H,J).rhs == -(0.5κa + 1.0im*ωa)*a_

# With coupling
ωa, ωb, κa, κb, g = (10.0,6.0,1.0,1.1,0.33)

H = ωa*(a'*a)⊗id + ωb*id⊗(b'*b) + g*(a'⊗b + a⊗b')
J = [sqrt(κa)*a⊗id,sqrt(κb)*id⊗b]

a_ = a⊗id
b_ = id⊗b
da = simplify_operators(1.0im*(H*a_ - a_*H))
da_sim = Qumulants.apply_comms(da)
da_qle = simplify_operators(1.0im*(H*a_ - a_*H) + sum(j'*a_*j - 0.5*(j'*j*a_ + a_*j'*j) for j=J))
da_qle_sim = Qumulants.apply_comms(da_qle)
@test da_qle_sim == heisenberg(a_,H,J).rhs == simplify_operators(-(0.5κa + 1.0im*ωa)*a_ - 1.0im*g*b_)

c = a'⊗b
dc = simplify_operators(1.0im*(H*c - c*H))
dc_sim = Qumulants.apply_comms(dc)
@test dc_sim == simplify_operators(1.0im*(ωa - ωb)*c + 1.0im*g*(b_'*b_ - a_'*a_))
dc_qle = simplify_operators(1.0im*(H*c - c*H) + sum(j'*c*j - 0.5*(j'*j*c + c*j'*j) for j=J))
dc_qle_sim = Qumulants.apply_comms(dc_qle)
@test dc_qle_sim == simplify_operators(1.0im*(ωa - ωb)*c + 1.0im*g*(b_'*b_ - a_'*a_) - 0.5*(κa + κb)*c)

end # testset
