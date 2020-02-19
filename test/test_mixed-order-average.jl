using Qumulants
using Test
using SymPy

@testset "mixed-order" begin

N = 2
# Parameters as SymPy variables
κ = symbols("κ", real=true, positive=true)
Δc = symbols("Δ_c", real=true)
Γ2 = [symbols("Γ^{($(k))}_2", real=true, positive=true) for k=1:N]
Γ3 = [symbols("Γ^{($(k))}_3", real=true, positive=true) for k=1:N]
Δ2 = [symbols("Δ^{($(k))}_2", real=true) for k=1:N]
Δ3 = [symbols("Δ^{($(k))}_3", real=true) for k=1:N]
Ω2 = [symbols("Ω^{($(k))}_2", real=true) for k=1:N]
Ω3 = [symbols("Ω^{($(k))}_3", real=true) for k=1:N]
g = [symbols("g^{($(k))}", real=true) for k=1:N]

# Operators
a = embed(Destroy(:a),1,N+1)
σ(i,j,k) = embed(Transition(Symbol("σ"*"^{($k)}"),i,j,3),k+1,N+1)



# Hamiltonian
H_atom = sum(-Δ2[k]*σ(2,2,k) - Δ3[k]*σ(3,3,k) + Ω2[k]*(σ(2,1,k) + σ(1,2,k)) + Ω3[k]*(σ(3,1,k) + σ(1,3,k)) for k=1:N)
H_cav = -Δc*a'*a + sum(g[k]*(a'*σ(1,2,k) + a*σ(2,1,k)) for k=1:N)
H = H_atom + H_cav;

# Jump operators
J = [sqrt(2κ)*a;[sqrt(Γ2[k])*σ(1,2,k) for k=1:N];[sqrt(Γ3[k])*σ(1,3,k) for k=1:N]]

# First-order operators
σ_ops = typeof(σ(1,2,1))[]
for k=1:N
    for i=1:3
        for j=i:3
            if !(i==j==1)
                push!(σ_ops, σ(i,j,k))
            end
        end
    end
end
ops_1st =[a; σ_ops];

# Second order products with the field
ops_2nd = AbstractOperator[a'*a]
for k=1:N
    for i=1:3
        for j=i:3
            if !(i==j==1)
                push!(ops_2nd, a'*σ(i,j,k))
                if i!=j
                    push!(ops_2nd, a*σ(i,j,k)) # RW terms ~ a*σ⁻
                end
            end
        end
    end
end

ops_tot = [ops_1st; ops_2nd];

# Derive second-order equations
eqs_2 = heisenberg(ops_tot,H,J);

orders = [2;ones(Int,N)]
eqs_avg = average(eqs_2,orders);

p = [κ;Δc;Γ2;Γ3;Δ2;Δ3;Ω2;Ω3;g]
meta_f = build_ode(eqs_avg,p;set_unknowns_zero=true)

# Test on which subspace stuff acts
tmp = σ(2,1,1)*σ(1,2,2)
tmp2 = eqs_2.rhs[13]
tmp3 = tmp2.args[6]
tmp4 = tmp3.args[1]
@test Qumulants.acts_on(tmp4) == Int[]
@test Qumulants.acts_on(tmp3) == [2,3]
@test Qumulants.acts_on(g[1]*a) == [1]
@test Qumulants.acts_on(tmp2) == [1,2,3]

end # testset
