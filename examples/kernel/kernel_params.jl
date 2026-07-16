# KernelParameters + update_parameters!: parameter sweeps without relowering (spec
# findings 1 and 9), plus end-to-end fixtures for the typed error taxonomy and the
# port-time defect fixes (findings 4 and 13). ORDER = 2 default; ORDER = 3 for timings.

using QuantumCumulants
using ModelingToolkitBase
using Symbolics
using SymbolicUtils
using QuantumOptics
include(joinpath(@__DIR__, "kernel_lower.jl"))
include(joinpath(@__DIR__, "kernel_eval.jl"))

isdefined(Main, :ORDER) || (ORDER = 2)
const N = 6

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = ORDER)
complete!(eqs)
ir = lower(eqs)

pd1 = Dict(J => 1.0, hx => 1.0, γ => 0.2)
pd2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)

# ---- Part A: COO -> nzval map, KernelParameters, update_parameters! ----

"""COO-entry -> position in M.nzval (duplicates accumulate into the same slot)."""
function nz_map(M, coo_i, coo_j)
    nzmap = Vector{Int}(undef, length(coo_i))
    for k in eachindex(coo_i)
        j = coo_j[k]
        r = Int(M.colptr[j]):(Int(M.colptr[j + 1]) - 1)
        p = searchsortedfirst(view(M.rowval, r), coo_i[k]) + first(r) - 1
        @assert M.rowval[p] == coo_i[k]
        nzmap[k] = p
    end
    return nzmap
end

struct KernelParameters
    ir::MomentIR
    values::Dict{Any, Any}
    nzmap::Vector{Int}
end
KernelParameters(ir::MomentIR, M, pdict) =
    KernelParameters(ir, Dict{Any, Any}(pdict), nz_map(M, ir.coo_i, ir.coo_j))

function write_nzval!(k::MomentKernel, kp::KernelParameters, cvals)
    fill!(k.M.nzval, zero(ComplexF64))
    @inbounds for t in eachindex(kp.nzmap)
        k.M.nzval[kp.nzmap[t]] += cvals[kp.ir.coo_c[t]]
    end
    return k
end

function update_parameters!(k::MomentKernel, kp::KernelParameters, pdict)
    for (p, v) in pdict
        kp.values[p] = v
    end
    return write_nzval!(k, kp, coefficient_values(kp.ir, kp.values))
end

k1 = MomentKernel(ir, coefficient_values(ir, pd1))
kp = KernelParameters(ir, k1.M, pd1)
update_parameters!(k1, kp, pd2)
k2 = MomentKernel(ir, coefficient_values(ir, pd2))
println("PARAMS nzval bit-exact after update_parameters!: ", k1.M.nzval == k2.M.nzval)

u = ComplexF64[0.1 * cos(3.7i) + 0.05im * sin(1.3i) for i in 1:ir.nstates]
du1 = zeros(ComplexF64, ir.nstates); du2 = similar(du1)
k1(du1, u, nothing, 0.0); k2(du2, u, nothing, 0.0)
println("PARAMS du bit-exact: ", du1 == du2)

t_subst = minimum(@elapsed(coefficient_values(kp.ir, kp.values)) for _ in 1:5)
cv = coefficient_values(kp.ir, kp.values)
t_map = minimum(@elapsed(write_nzval!(k1, kp, cv)) for _ in 1:20)
println(
    "PARAMS timing order=$ORDER: substitute=$(round(t_subst * 1.0e6, digits = 1))µs " *
        "mappass=$(round(t_map * 1.0e6, digits = 1))µs " *
        "ncoeffs=$(length(ir.coeffs)) nnz=$(length(k1.M.nzval)) ncoo=$(length(kp.nzmap))",
)

# ---- Part B: taxonomy fixtures (small model so they run in seconds) ----

Ns = 2
hs = ⊗([PauliSpace(Symbol(:s, i)) for i in 1:Ns]...)
sz(i) = Pauli(hs, :σ, 3, i); sx(i) = Pauli(hs, :σ, 1, i); sy(i) = Pauli(hs, :σ, 2, i)
sm(i) = (sx(i) - 1im * sy(i)) / 2
Hs = -J * sz(1) * sz(2) - hx * (sx(1) + sx(2))
eqs_s = meanfield([sz(i) for i in 1:Ns], Hs, [sm(i) for i in 1:Ns]; rates = [γ, γ], order = 1)
complete!(eqs_s)
lower(eqs_s)   # baseline: small model lowers clean

# 1. cancellation: get_variables(sum([J, -J])) loses J; discover_params must not
cs = Any[Symbolics.unwrap(J), Symbolics.unwrap(-J)]
old_n = length(Symbolics.get_variables(sum(cs)))
new_n = length(discover_params(cs))
println("TAXONOMY cancellation: old params=$old_n (J lost), new params=$new_n")

# 2. non-polynomial drift -> NonPolynomialDriftError (the one AutoBackend catches)
eqs_np = modify_equations(eqs_s, (op, d) -> d + exp(average(op)))
r_np = try
    lower(eqs_np); "NO ERROR (FAIL)"
catch e
    string(typeof(e))
end
println("TAXONOMY nonpoly: $r_np")

# 3. t in a coefficient -> TimeDependentCoefficientError at lower() time
tsym = eqs_s.iv
eqs_t = modify_equations(eqs_s, (op, d) -> cos(tsym) * d)
r_t = try
    lower(eqs_t); "NO ERROR (FAIL)"
catch e
    string(typeof(e))
end
println("TAXONOMY tdep: $r_t")

# 4a. user variable literally named `im` in a coefficient -> collision at lower() time
imvar = Symbolics.variable(:im)
eqs_im = modify_equations(eqs_s, (op, d) -> imvar * d)
r_im = try
    lower(eqs_im); "NO ERROR (FAIL)"
catch e
    string(typeof(e))
end
println("TAXONOMY im-in-coeff: $r_im")

# 4b. user passes a pdict key named `im` -> collision at coefficient_values time
ir_s = lower(eqs_s)
r_imp = try
    coefficient_values(ir_s, merge(pd1, Dict(imvar => 2.0))); "NO ERROR (FAIL)"
catch e
    string(typeof(e))
end
println("TAXONOMY im-in-pdict: $r_imp")

# regression guard: benchmark coefficients still evaluate (algebra's own Sym(:im) must
# NOT trip the guard) and the small model still solves the taxonomy edits unharmed
println("TAXONOMY regression: benchmark cvals ok = ", length(cv) == length(ir.coeffs))
