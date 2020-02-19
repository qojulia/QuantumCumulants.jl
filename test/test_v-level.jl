using Qumulants
using SymPy
using Test
import OrdinaryDiffEq
ODE = OrdinaryDiffEq

@testset "v-level" begin

# Parameters as SymPy variables
Γ2, Γ3, κ = symbols(" Γ_2 Γ_3 κ", real=true, positive=true)
Δ2, Δ3, Ω2, Ω3, Δc, g = symbols("Δ_2 Δ_3 Ω_2 Ω_3 Δ_c g", real=true)

a = Destroy(:a) ⊗ Identity()
σ(i,j) = Identity() ⊗ Transition(:σ,i,j,3)

for i=2:3
    @test average(σ(i,i)).__pyobject__.is_real
end

H_atom = -Δ2*σ(2,2) - Δ3*σ(3,3) + Ω2*(σ(2,1) + σ(1,2)) + Ω3*(σ(3,1) + σ(1,3))
H_cav = -Δc*a'*a + g*(a'*σ(1,2) + a*σ(2,1))
H = H_atom + H_cav

J = [sqrt(2κ)*a,sqrt(Γ2)*σ(1,2),sqrt(Γ3)*σ(1,3)]

# Operators for mean-field eqs
σ_ops = typeof(σ(1,2))[]
for i=1:3
    for j=i:3
        if !(i==j==1)
            push!(σ_ops, σ(i,j))
        end
    end
end
ops_mf =[a; σ_ops]

# Derive eqs
op_eqs = heisenberg(ops_mf,H,J)

eqs_mf = average(op_eqs,1)

p = [Γ2, Γ3, κ, Δ2, Δ3, Ω2, Ω3, Δc, g]
meta_mf = build_ode(eqs_mf,p)

# Second order equations
ops_2nd = [a'*a;ops_mf]
for i=1:3
    for j=i:3
        if !(i==j==1)
            push!(ops_2nd, a'*σ(i,j))
            if i!=j
                push!(ops_2nd, a*σ(i,j))
            end
        end
    end
end

ops2_eqs = heisenberg(ops_2nd,H,J)

eqs_2nd = average(ops2_eqs,2)

meta_f = build_ode(eqs_2nd,p;set_unknowns_zero=true)

end # testset
