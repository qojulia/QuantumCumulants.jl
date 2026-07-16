# Kernel cache (issue #294, spec "Kernel cache design (v1)").
#
# Key = SHA-256 of a canonical TEXT serialization (version stamps + equations in state
# order). The digest is ONLY an index: on a hit the stored text is compared
# byte-for-byte, so a false hit degrades to a miss. Payload = plain data only (Int32
# tables, parameter names, and the pooled coefficient evaluator as PRINTED JULIA SOURCE,
# turned into an RGF on load), never serialized symbolic objects. A loaded kernel is
# therefore fully sweepable. Corrupt/unreadable files warn and fall back to lowering.
#
# Storage lives in the JLD2 package extension (`QuantumCumulantsJLD2Ext`), reached
# through `_CACHE_BACKEND`; without `using JLD2` the cache keyword is a typed error.

const IR_FORMAT_VERSION = 1

_pkgver(m) = try
    string(pkgversion(m))
catch
    "unknown"
end

"""
Canonical text of a completed equation set: version stamps, then every `state := drift`
line in `keys(eqs.graph.nodes)` order. Parameters are NOT part of the text (they live in
the payload), so the digest is computable before lowering.
"""
function canonical_text(eqs)
    io = IOBuffer()
    println(io, "ir_format=", IR_FORMAT_VERSION)
    println(io, "julia=", VERSION)
    println(io, "QuantumCumulants=", _pkgver(@__MODULE__))
    println(io, "SecondQuantizedAlgebra=", _pkgver(SQA))
    println(io, "SymbolicUtils=", _pkgver(SymbolicUtils))
    println(io, "Symbolics=", _pkgver(Symbolics))
    for (k, nd) in pairs(eqs.graph.nodes)
        println(io, string(k), " := ", string(nd.drift))
    end
    return String(take!(io))
end

"""
Pooled coefficient evaluator as printed Julia source: f(pvec) -> coefficient vector.
The algebra's `Sym(:im)` prints as the literal `im`, which resolves to `Base.im` inside
the RGF module (verified numerically in the prototype).
"""
coeff_source(ir) = string(
    Symbolics.build_function(
        map(Symbolics.wrap, ir.coeffs), map(Symbolics.wrap, ir.params);
        expression = Val{true},
    )[1],
)

# The JLD2 extension installs `(store! = ..., load = ...)` here in its `__init__`.
const _CACHE_BACKEND = Ref{Any}(nothing)

function _cache_backend()
    b = _CACHE_BACKEND[]
    b === nothing && throw(
        ArgumentError("KernelBackend(cache = path) requires JLD2; run `using JLD2` first.")
    )
    return b
end

function _cache_store!(path, text, digest, ir)
    payload = Dict{String, Any}(
        "nstates" => ir.nstates,
        "parent" => ir.parent,
        "leaf" => ir.leaf,
        "coo_i" => ir.coo_i,
        "coo_j" => ir.coo_j,
        "coo_c" => ir.coo_c,
        "params" => String[string(p) for p in ir.params],
        "coeff_src" => coeff_source(ir),
    )
    _cache_backend().store!(path, text, digest, payload)
    return nothing
end

"""Returns the loaded entry or `nothing` (miss / failed verify / unreadable file)."""
function _cache_load(path, text, digest)
    raw = _cache_backend().load(path, text, digest)
    raw === nothing && return nothing
    return (
        nstates = raw["nstates"]::Int,
        parent = raw["parent"]::Vector{Int32},
        leaf = raw["leaf"]::Vector{Int32},
        coo_i = raw["coo_i"]::Vector{Int32},
        coo_j = raw["coo_j"]::Vector{Int32},
        coo_c = raw["coo_c"]::Vector{Int32},
        params = raw["params"]::Vector{String},
        rgf = @RuntimeGeneratedFunction(Meta.parse(raw["coeff_src"]::String)),
    )
end
