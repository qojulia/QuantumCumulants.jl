# Evaluation: MomentIR -> in-place RHS callable (issue #294, M·v design).
#
# The RHS is two data passes, no generated code:
#   1. update the distinct-monomial vector v via prefix chains (one fused gather-multiply
#      per monomial; parents have smaller ids, so ascending order is a valid schedule),
#   2. du = M * v (one complex SpMV).
#
# Note: the kernel carries its workspace `v`, so one kernel instance must not be called
# concurrently from multiple threads (same caveat as any cache-carrying RHS closure).

struct MomentKernel
    M::SparseMatrixCSC{ComplexF64, Int32}
    parent::Vector{Int32}
    leaf::Vector{Int32}
    v::Vector{ComplexF64}
end

# `ir` untyped so this constructor also serves cache hits, which carry the plain tables
# (M, parent, leaf) without ever constructing a MomentIR.
function MomentKernel(ir, cvals::Vector{ComplexF64})
    v = zeros(ComplexF64, length(ir.parent))
    v[1] = one(ComplexF64)
    return MomentKernel(assemble(ir, cvals), ir.parent, ir.leaf, v)
end

"""Refresh the distinct-monomial vector in place (shared by the RHS and the Jacobian)."""
function update_v!(v, parent, leaf, u)
    @inbounds for m in 2:length(v)
        j = leaf[m]
        x = j > 0 ? u[j] : conj(u[-j])
        v[m] = v[parent[m]] * x
    end
    return v
end

function (k::MomentKernel)(du, u, p, t)
    update_v!(k.v, k.parent, k.leaf, u)
    mul!(du, k.M, k.v)
    return nothing
end

# ---- parameter sweeps without relowering ----------------------------------------------

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

"""
The parameter payload of a kernel `ODEProblem` (`prob.p`). Carries the discovered
parameter occurrences, the current values, and everything needed to rewrite `M.nzval`
in place on a parameter update. `evalcoeffs` abstracts how pooled coefficients are
evaluated: fresh lowerings substitute into the symbolic coefficients, cache-loaded
kernels call the stored coefficient evaluator.
"""
struct KernelParameters{F}
    params::Vector{Any}
    values::Dict{Any, Any}
    evalcoeffs::F                 # values dict -> Vector{ComplexF64} of pooled coefficients
    coo_c::Vector{Int32}
    nzmap::Vector{Int}
end

function KernelParameters(ir::MomentIR, M::SparseMatrixCSC, values::Dict)
    evalcoeffs = vals -> coefficient_values(ir, vals)
    return KernelParameters(
        ir.params, Dict{Any, Any}(values), evalcoeffs, ir.coo_c,
        nz_map(M, ir.coo_i, ir.coo_j),
    )
end

function write_nzval!(k::MomentKernel, kp::KernelParameters, cvals)
    fill!(k.M.nzval, zero(ComplexF64))
    @inbounds for t in eachindex(kp.nzmap)
        k.M.nzval[kp.nzmap[t]] += cvals[kp.coo_c[t]]
    end
    return k
end
