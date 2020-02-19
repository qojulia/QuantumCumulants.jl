using Qumulants
using Test
using BenchmarkTools
using PyPlot
using InteractiveUtils

# Test single mode
a = Destroy(:a)

ω, κ = (10.0, 0.1)
H = ω*a'*a
J = [sqrt(κ)*a]
Jdagger = adjoint.(J)

da = 1.0im*(H*a - a*H)
@test heisenberg(a,H) == da == -1.0im*ω*a
da_qle = 1.0im*(H*a - a*H) + sum(j'*a*j - 0.5*(j'*j*a + a*j'*j) for j=J)
@test heisenberg(a,H,J) == da_qle == -(1.0im*ω + 0.5κ)*a

bench = @benchmark heisenberg($a,$H,$J;Jdagger=$Jdagger)

# Test composite system
b = Destroy(:b)
id = Identity()

ωa, ωb, κa, κb, g = (10.0,6.0,1.0,1.1,0.5)

H = ωa*(a'*a)⊗id + ωb*id⊗(b'*b) + g*(a'⊗b + a⊗b')
c = a⊗b
J = [sqrt(κa)*a⊗id,sqrt(κb)*id⊗b]
Jdagger = adjoint.(J)

a_ = a⊗id
da = 1.0im*(H*a_ - a_*H)

bench2 = @benchmark heisenberg($a_,$H,$J;Jdagger=$Jdagger)


function bench_time(N)
    labels = [Symbol("a$i") for i=1:N]
    ops = [embed(Destroy(labels[i]), (1,N), i) for i=1:N]
    ω = ones(N)
    κ = ω./2
    g = ones(N,N)
    for i=1:N
        g[i,i] = 0.0
    end

    H = sum(ω[i]*ops[i]'*ops[i] for i=1:N) + sum(g[i,j]*ops[i]'*ops[j] for i=1:N, j=1:N)
    J = ops
    Jdagger = adjoint.(J)

    a = ops[1]
    da = heisenberg(a,H,J;Jdagger=Jdagger,rates=κ)
    t = @belapsed heisenberg($a,$H,$J;Jdagger=$Jdagger,rates=$κ) samples=3 evals=5
    return t
end


N_ls = [2:6;]
times = Float64[]
for N=N_ls
    println("Benchmarking $N")
    push!(times, bench_time(N))
end

semilogy(N_ls, times, "o", label="Single Eq. Time")
# plot(N_ls, N_ls .* times, "--", label="Extrapolated to N Equations")
