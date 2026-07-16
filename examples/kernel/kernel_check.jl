# Correctness of the M·v kernel against a codegen-free reference: numerically substitute
# random (u, p) into the symbolic drifts and compare du equation by equation. Safe to run
# in one session (no compile timing here).

using QuantumCumulants
using Symbolics
using SymbolicUtils
using Random
include(joinpath(@__DIR__, "kernel_lower.jl"))
include(joinpath(@__DIR__, "kernel_eval.jl"))

function check_du(eqs, pdict; nsample = 40, npoints = 3, seed = 1)
    g = eqs.graph
    ir = lower(eqs)
    k = MomentKernel(ir, coefficient_values(ir, pdict))
    vars, idx = statevars(g)
    drifts = [Symbolics.unwrap(nd.drift) for nd in values(g.nodes)]
    n = length(drifts)
    rng = MersenneTwister(seed)
    sel = n <= nsample ? collect(1:n) : sort!(randperm(rng, n)[1:nsample])
    pd = Dict(Symbolics.unwrap(kk) => vv for (kk, vv) in pdict)
    maxrel = 0.0
    for _ in 1:npoints
        u = randn(rng, ComplexF64, n)
        du = similar(u)
        k(du, u, nothing, 0.0)
        subs = Dict{Any, Any}(pd)
        for v in vars
            j = idx[v]
            subs[v] = j > 0 ? u[j] : conj(u[-j])
        end
        # the symbolic imaginary unit appears in the drifts as well (see kernel_lower.jl)
        for i in sel, v in Symbolics.get_variables(drifts[i])
            uu = SymbolicUtils.unwrap(v)
            SymbolicUtils.issym(uu) && SymbolicUtils.nameof(uu) === :im && (subs[uu] = im)
        end
        for i in sel
            ref = ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(drifts[i], subs)))
            maxrel = max(maxrel, abs(du[i] - ref) / max(abs(ref), 1e-12))
        end
    end
    return maxrel
end

# Pauli chain, order 2 (all equations checked)
N = 3
h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = 2)
complete!(eqs)
e1 = check_du(eqs, Dict(J => 1.0, hx => 1.0, γ => 0.2))
println("pauli N=3 order=2 (36 eqs, all):        max rel err = ", e1)

# Kerr cavity, conj-folded (partner resolution, powers, complex coefficients)
hc = FockSpace(:cavity)
a = Destroy(hc, :a)
@variables Δ Ω κ U
Hk = Δ * a' * a + Ω * (a + a') + 0.5 * U * (a' * a' * a * a)
eqs_k = meanfield([a], Hk, [a]; rates = [κ], order = 2)
complete!(eqs_k; get_adjoints = false)
e2 = check_du(eqs_k, Dict(Δ => 0.5, Ω => 1.3, κ => 0.7, U => 0.4))
println("kerr folded order=2 (3 eqs, all):        max rel err = ", e2)

# Benchmark system, sampled
N6 = 6
h6 = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N6]...)
sx(i) = Pauli(h6, :σ, 1, i); sy(i) = Pauli(h6, :σ, 2, i); sz(i) = Pauli(h6, :σ, 3, i)
sm(i) = (sx(i) - 1im * sy(i)) / 2
H6 = -J * sum(sz(i) * sz(i + 1) for i in 1:(N6 - 1)) - hx * sum(sx(i) for i in 1:N6)
eqs6 = meanfield([sz(i) for i in 1:N6], H6, [sm(i) for i in 1:N6]; rates = [γ for i in 1:N6], order = 3)
complete!(eqs6)
e3 = check_du(eqs6, Dict(J => 1.0, hx => 1.0, γ => 0.2))
println("ising N=6 order=3 (693 eqs, 40 sampled): max rel err = ", e3)

ok = max(e1, e2, e3) < 1e-10
println(ok ? "PASS" : "FAIL")
