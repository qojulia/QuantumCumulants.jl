# Safe kernel cache prototype (spec "Kernel cache design (v1)", finding 2).
#
# Key = SHA-256 of a canonical TEXT serialization (equations in state order + params +
# version stamps). The digest is ONLY an index: on a hit the stored text is compared
# byte-for-byte, so a false hit degrades to a miss. Payload = plain data only (Int32
# tables, parameter names, and the pooled coefficient evaluator as PRINTED JULIA SOURCE,
# turned into an RGF on load), never serialized symbolic objects. A loaded kernel is
# therefore fully sweepable. Corrupt/unreadable files warn and fall back to lowering.
#
#   ORDER = 2 (default) | 3;  MODE = :all (default, full test) | :digest (print digest only,
#   for the two-process cross-session stability check)

using QuantumCumulants
using ModelingToolkitBase
using Symbolics
using SymbolicUtils
using QuantumOptics
using SHA
using JLD2
using SparseArrays
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
include(joinpath(@__DIR__, "kernel_lower.jl"))
include(joinpath(@__DIR__, "kernel_eval.jl"))

const IR_FORMAT_VERSION = 1

isdefined(Main, :ORDER) || (ORDER = 2)
isdefined(Main, :MODE) || (MODE = :all)
const N = 6

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = ORDER)
complete!(eqs)
ir = lower(eqs)

_pkgver(m) = try string(pkgversion(m)) catch; "unknown" end

function canonical_text(eqs, ir)
    io = IOBuffer()
    println(io, "ir_format=", IR_FORMAT_VERSION)
    println(io, "julia=", VERSION)
    println(io, "QuantumCumulants=", _pkgver(QuantumCumulants))
    println(io, "SymbolicUtils=", _pkgver(SymbolicUtils))
    println(io, "Symbolics=", _pkgver(Symbolics))
    println(io, "params=", join(string.(ir.params), ","))
    for (k, nd) in pairs(eqs.graph.nodes)
        println(io, string(k), " := ", string(nd.drift))
    end
    return String(take!(io))
end

"""Pooled coefficient evaluator as printed Julia source: f(pvec) -> coefficient vector.
The algebra's `Sym(:im)` prints as the literal `im`, which resolves to `Base.im` inside
the RGF module (verified numerically below)."""
coeff_source(ir) = string(
    Symbolics.build_function(
        map(Symbolics.wrap, ir.coeffs), map(Symbolics.wrap, ir.params);
        expression = Val{true},
    )[1],
)

function cache_store!(path, text, digest, ir)
    jldopen(path, "a+") do f
        haskey(f, digest) && return nothing
        f["$digest/text"] = text
        f["$digest/nstates"] = ir.nstates
        f["$digest/parent"] = ir.parent
        f["$digest/leaf"] = ir.leaf
        f["$digest/coo_i"] = ir.coo_i
        f["$digest/coo_j"] = ir.coo_j
        f["$digest/coo_c"] = ir.coo_c
        f["$digest/params"] = string.(ir.params)
        f["$digest/coeff_src"] = coeff_source(ir)
        return nothing
    end
    return nothing
end

"""Returns the loaded entry or `nothing` (miss / failed verify / unreadable file)."""
function cache_load(path, text, digest)
    isfile(path) || return nothing
    return try
        jldopen(path, "r") do f
            haskey(f, digest) || return nothing
            f["$digest/text"] == text || return nothing   # byte-exact verify
            (
                nstates = f["$digest/nstates"],
                parent = f["$digest/parent"], leaf = f["$digest/leaf"],
                coo_i = f["$digest/coo_i"], coo_j = f["$digest/coo_j"],
                coo_c = f["$digest/coo_c"], params = f["$digest/params"],
                rgf = @RuntimeGeneratedFunction(Meta.parse(f["$digest/coeff_src"])),
            )
        end
    catch err
        @warn "kernel cache unreadable; lowering fresh" path err
        nothing
    end
end

function loaded_kernel(e, pdict)
    byname = Dict(string(Symbolics.unwrap(k)) => v for (k, v) in pdict)
    cvals = ComplexF64.(e.rgf([byname[n] for n in e.params]))
    M = sparse(e.coo_i, e.coo_j, cvals[e.coo_c], e.nstates, length(e.parent), +)
    v = zeros(ComplexF64, length(e.parent)); v[1] = one(ComplexF64)
    return MomentKernel(M, e.parent, e.leaf, v)
end

text = canonical_text(eqs, ir)
digest = bytes2hex(sha256(text))
println("DIGEST order=$ORDER $digest")
MODE === :digest && exit(0)

path = joinpath(@__DIR__, "kernel_cache_test.jld2")
isfile(path) && rm(path)
cache_store!(path, text, digest, ir)
println("CACHE file size = $(filesize(path)) bytes")

pd1 = Dict(J => 1.0, hx => 1.0, γ => 0.2)
pd2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)

# 1. round-trip at pd1
t_load = @elapsed e = cache_load(path, text, digest)
@assert e !== nothing
k_f = MomentKernel(ir, coefficient_values(ir, pd1))
k_l = loaded_kernel(e, pd1)
u = ComplexF64[0.1 * cos(3.7i) + 0.05im * sin(1.3i) for i in 1:ir.nstates]
du_f = zeros(ComplexF64, ir.nstates); du_l = similar(du_f)
k_f(du_f, u, nothing, 0.0); k_l(du_l, u, nothing, 0.0)
dnz = k_f.M.nzval == k_l.M.nzval
println("CACHE roundtrip: nzval bit-exact=$(dnz) du bit-exact=$(du_f == du_l)",
    dnz ? "" : " (max rel dev $(maximum(abs.(k_f.M.nzval .- k_l.M.nzval) ./ abs.(k_f.M.nzval))))")

# 2. sweep after load: loaded evaluator at pd2 vs fresh substitute path at pd2
k_f2 = MomentKernel(ir, coefficient_values(ir, pd2))
k_l2 = loaded_kernel(e, pd2)
println("CACHE sweep-after-load: nzval bit-exact=", k_f2.M.nzval == k_l2.M.nzval)

# 3. simulated collision: right digest, tampered stored text -> must be a miss
path2 = joinpath(@__DIR__, "kernel_cache_tamper.jld2")
isfile(path2) && rm(path2)
jldopen(path2, "w") do f
    f["$digest/text"] = text * "TAMPERED"
end
println("CACHE tampered-text treated as miss: ", cache_load(path2, text, digest) === nothing)
rm(path2)

# 4. corrupt file -> warn + nothing (caller lowers fresh)
path3 = joinpath(@__DIR__, "kernel_cache_corrupt.jld2")
write(path3, rand(UInt8, 512))
println("CACHE corrupt-file degrades to miss: ", cache_load(path3, text, digest) === nothing)
rm(path3)

# 5. load-path timing (open + verify + parse/RGF + eval + assemble)
t_full = @elapsed begin
    e2 = cache_load(path, text, digest)
    loaded_kernel(e2, pd1)
end
println("CACHE load-path: first_load=$(round(t_load, digits = 3))s full_reload=$(round(t_full, digits = 3))s")
