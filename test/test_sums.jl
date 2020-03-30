using Qumulants
using Test
using SymPy: symbols
using MacroTools
using OrdinaryDiffEq

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
H_TC = Sum(Δ[i]*σ(:e,:e,i),i) + Sum(g[i]*(a'*σ(:g,:e,i) + a*σ(:e,:g,i)),i)
J = [a,σ(:g,:e,i),σ(:e,:g,i)]

ops = [(a'*a),a'*σ(:g,:e,j),σ(:e,:e,j),σ(:e,:g,k)*σ(:g,:e,j)]
he = heisenberg(ops,H_TC,J;rates=[2κ,γ[i],ν[i]])

he_avg = average(he,2)

# tmp = Qumulants.order_indexed(he_avg.rhs,ops,2)

n = symbols("n", integer=true)
p = (g,Δ,κ,γ,ν,n-1)
meta_f = build_ode(he_avg,p;set_unknowns_zero=true)
f = Meta.eval(meta_f)


# Implementation by hand
function f2(du,u,t,p)
    # Assign stuff
    n = u[1,1]
    adσ = @view u[1,2:end]
    dadσ = @view du[1,2:end]
    σpσ = @view u[2:end,2:end]
    dσpσ = @view du[2:end,2:end]

    (g,Δ,κ,γ,ν,N) = p0

    # Photon number
    du[1,1] = -2κ*n

    # σ_i⁺σ_j
    for j=1:N

        # Add terms to photon number
        du[1,1] += 2.0*g[j]*imag(adσ[j])

        # adσ
        dadσ[j] = -(1.0im*Δ[j] + κ + (ν[j] + γ[j])/2.0)*adσ[j] + 1.0im*g[j]*(2σpσ[j,j] - 1.0)*n
        for m=1:N
            dadσ[j] += 1.0im*g[m]*σpσ[m,j]
        end

        # Diagonal σpσm elements
        dσpσ[j,j] = -(γ[j] + ν[j])*σpσ[j,j] + ν[j] - 2.0*g[j]*imag(adσ[j])

        # Off-diagonal elements
        for k=j+1:N
            dσpσ[j,k] = ( (1.0im*(Δ[j] - Δ[k]) - (ν[j] + ν[k] + γ[j] + γ[k])/2.0)*σpσ[j,k] +
                            -1.0im*g[j]*(2σpσ[j,j] - 1.0)*adσ[k] + 1.0im*g[k]*(2σpσ[k,k] - 1.0)*conj(adσ[j])
                        )
        end

        # Lower triangle of σpσ is trivial
        for k=1:j-1
            dσpσ[j,k] = conj(dσpσ[k,j])
        end

        for k=1:N
            du[k+1,1] = conj(dadσ[k])
        end

    end
end

meta_tmp = :( sum((-0.5 * Int(i != j) * Int(j != i) * Int(k != j) * (p[5])[i] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((0.5 * Int(i != j) * Int(j != i) * Int(k != i) * Int(k != j) * (p[5])[i] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-0.5 * Int(i != j) ^ 2 * Int(j != i) * Int(k != j) * (p[4])[i] * u[i + p[6] + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((0.5 * Int(i != j) ^ 2 * Int(j != i) * Int(k != j) * (p[5])[i] * u[i + p[6] + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-1.0 * im * Int(i != j) ^ 2 * Int(j != i) * Int(k != j) * (p[1])[i] * u[i + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-1.0 * im * Int(i != j) ^ 2 * Int(j != i) * Int(k != j) * (p[1])[i] * u[j + 1] * u[k + p[6] * (i - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-1.0 * im * Int(i != j) ^ 2 * Int(j != i) * Int(k != j) * (p[2])[i] * u[i + p[6] + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-1.0 * im * Int(i != j) ^ 2 * Int(j != i) * adjoint(Int(k != j)) * adjoint(u[k + 1]) * (p[1])[i] * u[i + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((0.5 * Int(i != j) ^ 2 * Int(j != i) * Int(k != i) * Int(k != j) * (p[4])[i] * u[i + p[6] + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-0.5 * Int(i != j) ^ 2 * Int(j != i) * Int(k != i) * Int(k != j) * (p[5])[i] * u[i + p[6] + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((-1.0 * im * Int(i != j) * Int(j != i) * Int(k != j) * adjoint(Int(i != j)) * adjoint(u[i + 1]) * (p[1])[i] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((1.0 * im * Int(i != j) ^ 2 * Int(j != i) * Int(k != i) * Int(k != j) * (p[1])[i] * u[i + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((1.0 * im * Int(i != j) ^ 2 * Int(j != i) * Int(k != i) * Int(k != j) * (p[1])[i] * u[j + 1] * u[k + p[6] * (i - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((1.0 * im * Int(i != j) ^ 2 * Int(j != i) * Int(k != i) * Int(k != j) * (p[2])[i] * u[i + p[6] + 1] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((1.0 * im * Int(i != j) ^ 2 * Int(j != i) * adjoint(Int(k != i)) * adjoint(Int(k != j)) * adjoint(u[k + 1]) * (p[1])[i] * u[i + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])) + sum((1.0 * im * Int(i != j) * Int(j != i) * Int(k != i) * Int(k != j) * adjoint(Int(i != j)) * adjoint(u[i + 1]) * (p[1])[i] * u[k + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6])))
using MacroTools
MacroTools.postwalk(x->MacroTools.@capture(x, Sumsym_(arg_ for i_=l_:u_)) ? 0 : x, meta_tmp)

function f3(du, u, p, t)
          begin
              begin
                  du[1] = -2.0 * p[3] * u[1] + sum((1.0 * im * adjoint(u[i + 1]) * (p[1])[i] for i = 1:p[6])) + sum((-1.0 * im * (p[1])[i] * u[i + 1] for i = 1:p[6]))
                  begin
                      for j = 1:p[6]
                          du[j + 1] = ((((((-1.0 * p[3] * u[j + 1] + 2.0 * im * u[1] * (p[1])[j] * u[j + p[6] + 1]) - 1.0 * im * u[1] * (p[1])[j]) + 1.0 * im * (p[1])[j] * u[j + p[6] + 1]) - 1.0 * im * (p[2])[j] * u[j + 1]) - 0.5 * (p[4])[j] * u[j + 1]) - 0.5 * (p[5])[j] * u[j + 1]) + sum((1.0 * im * Int(i != j) ^ 2 * Int(j != i) * (p[1])[i] * u[i + p[6] * (j - 1) + 2 * p[6] + 1] for i = 1:p[6]))
                          du[j + p[6] + 1] = (((-1.0 * im * adjoint(u[j + 1]) * (p[1])[j] + 1.0 * im * (p[1])[j] * u[j + 1]) - 1.0 * (p[4])[j] * u[j + p[6] + 1]) - 1.0 * (p[5])[j] * u[j + p[6] + 1]) + 1.0 * (p[5])[j]
                          for k = 1:p[6]
                              k != j && (du[k + p[6] * (j - 1) + 2 * p[6] + 1] = ((((((((-1.0 * im * Int(k != j) ^ 2 * (p[1])[k] * u[j + 1] * u[k + p[6] + 1] + 1.0 * im * Int(k != j) ^ 2 * (p[2])[k] * u[k + p[6] * (j - 1) + 2 * p[6] + 1]) - 0.5 * Int(k != j) ^ 2 * (p[4])[k] * u[k + p[6] * (j - 1) + 2 * p[6] + 1]) + 1.0 * im * Int(k != j) * (p[1])[k] * u[j + 1]) - 1.0 * im * Int(k != j) * (p[2])[j] * u[k + p[6] * (j - 1) + 2 * p[6] + 1]) - 0.5 * Int(k != j) * (p[4])[j] * u[k + p[6] * (j - 1) + 2 * p[6] + 1]) - 0.5 * Int(k != j) * (p[5])[j] * u[k + p[6] * (j - 1) + 2 * p[6] + 1]) + 2.0 * im * adjoint(Int(k != j)) * adjoint(u[k + 1]) * (p[1])[j] * u[j + p[6] + 1]) - 1.0 * im * adjoint(Int(k != j)) * adjoint(u[k + 1]) * (p[1])[j]) )
                          end
                      end
                  end
              end
          end
          return nothing
end


# Test numerical result
N = 2

κn = 1
Δn = zeros(N)
gn = ones(N)
γn = zeros(N)
νn = 3 .* ones(N)
p0 = (gn,Δn,κn,γn,νn,N)
tmax = 10.0

u0 = zeros(ComplexF64,(N+1)^2)
u2 = zeros(ComplexF64, 1+N,1+N)

prob = ODEProblem(f,u0,(0.0,tmax),p0)
prob2 = ODEProblem(f2,u2,(0.0,tmax),p0)
prob3 = ODEProblem(f3,u0,(0.0,tmax),p0)

sol = solve(prob,Tsit5())
sol2 = solve(prob2,Tsit5())
sol3 = solve(prob3,Tsit5())
