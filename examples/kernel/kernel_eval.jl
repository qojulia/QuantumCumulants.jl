# Evaluation: MomentIR -> in-place RHS callable + ODEProblem (issue #294, M·v design).
#
# The RHS is two data passes, no generated code:
#   1. update the distinct-monomial vector v via prefix chains (one fused gather-multiply
#      per monomial; parents have smaller ids, so ascending order is a valid schedule),
#   2. du = M * v (one complex SpMV).
#
# Note: the kernel carries its workspace `v`, so one kernel instance must not be called
# concurrently from multiple threads (same caveat as any cache-carrying RHS closure).

using SparseArrays
using LinearAlgebra: mul!

struct MomentKernel
    M::SparseMatrixCSC{ComplexF64, Int32}
    parent::Vector{Int32}
    leaf::Vector{Int32}
    v::Vector{ComplexF64}
end

# `ir` untyped so this file is loadable without kernel_lower.jl (cross-session warm starts
# deserialize (M, parent, leaf) directly and never construct a MomentIR)
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

"""
    moment_kernel(eqs, pdict) -> MomentKernel

Lower a completed `MeanfieldEquations` and materialize the kernel for one parameter
assignment. Parameter sweeps: `coefficient_values(ir, pdict)` + `MomentKernel(ir, c)`
(or write `kernel.M.nzval` in place through `assemble`'s COO order).
"""
function moment_kernel(eqs, pdict)
    ir = lower(eqs)
    return MomentKernel(ir, coefficient_values(ir, pdict))
end
