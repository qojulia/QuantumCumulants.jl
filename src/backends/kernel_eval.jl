# Evaluation: MomentIR -> in-place RHS callable (issue #294, M·v design).
#
# The RHS is two data passes, no generated code:
#   1. update the distinct-monomial vector v via prefix chains (one fused gather-multiply
#      per monomial; parents have smaller ids, so ascending order is a valid schedule),
#   2. du = M * v (one complex SpMV).
#
# Note: the kernel carries its workspace `v`, so one kernel instance must not be called
# concurrently from multiple threads (same caveat as any cache-carrying RHS closure).

# Equation count above which threading the RHS pays off end to end. Both passes thread
# through Polyester's persistent task pool (`@batch`), whose per-call dispatch cost is ~1 µs
# (vs ~7 µs for `Threads.@threads`), so the full RHS wins from ~65 equations up (measured
# 1.1x at ~65 eqs, 2.4-3.1x by a few hundred, 12 threads); below that the pool overhead
# loses on a sub-microsecond RHS. Only consulted for `parallel = :auto`.
const KERNEL_PARALLEL_MIN = 64
_resolve_kernel_parallel(flag, neq) =
    flag === :auto ? (Threads.nthreads() > 1 && neq >= KERNEL_PARALLEL_MIN) : flag

# The RHS stores Mᵀ (the coefficient matrix transposed, i.e. M in CSR): `du = M * v` then
# becomes a row-wise gather where each `du[i]` sums one column of `Mt` independently. That
# layout is both cache-friendlier than a CSC `mul!` (sequential `du` writes) and trivially
# parallel (no write conflicts). `fac`/`fac_ptr` are the flat factor chains used by the
# threaded monomial update (see `update_v_flat!`); they are derived from `parent`/`leaf` and
# add no cache format. The kernel carries its workspace `v`, so one instance must not be
# called concurrently from multiple threads (same caveat as any cache-carrying RHS).
struct MomentKernel
    Mt::SparseMatrixCSC{ComplexF64, Int32}   # transpose of the coefficient matrix M
    parent::Vector{Int32}
    leaf::Vector{Int32}
    fac::Vector{Int32}                       # flat factor chains (threaded update path)
    fac_ptr::Vector{Int32}                   # CSR offsets into `fac`, one range per monomial
    v::Vector{ComplexF64}
    parallel::Bool                           # thread the RHS (monomial update + SpMV)
end

# `ir` untyped so this constructor also serves cache hits, which carry the plain tables
# (parent, leaf) without ever constructing a MomentIR.
function MomentKernel(ir, cvals::Vector{ComplexF64}; parallel::Bool = false)
    v = zeros(ComplexF64, length(ir.parent))
    v[1] = one(ComplexF64)
    fac, fac_ptr = build_flat(ir.parent, ir.leaf)
    return MomentKernel(assemble(ir, cvals), ir.parent, ir.leaf, fac, fac_ptr, v, parallel)
end

# Convenience constructor over prebuilt tables (cache load / copy): derives `fac`/`fac_ptr`.
function MomentKernel(
        Mt::SparseMatrixCSC, parent::Vector{Int32}, leaf::Vector{Int32},
        v::Vector{ComplexF64}, parallel::Bool,
    )
    fac, fac_ptr = build_flat(parent, leaf)
    return MomentKernel(Mt, parent, leaf, fac, fac_ptr, v, parallel)
end

"""Refresh the distinct-monomial vector in place via the prefix chains (serial; also used by
the Jacobian). Each monomial is `(parent monomial) * (one state factor)`, and parents have
smaller ids, so ascending order is a valid schedule."""
function update_v!(v, parent, leaf, u)
    @inbounds for m in 2:length(v)
        j = leaf[m]
        x = j > 0 ? u[j] : conj(u[-j])
        v[m] = v[parent[m]] * x
    end
    return v
end

"""Flat factor chains: `fac[fac_ptr[m]:fac_ptr[m+1]-1]` is the full (signed) factor list of
monomial `m`, in the same order the prefix chain would append them. Storing the whole chain
makes each monomial an INDEPENDENT product, so `update_v_flat!` threads without the prefix
form's cross-iteration dependency; the shared order keeps it bit-identical to `update_v!`."""
function build_flat(parent::Vector{Int32}, leaf::Vector{Int32})
    n = length(parent)
    depth = zeros(Int32, n)
    @inbounds for m in 2:n
        depth[m] = depth[parent[m]] + Int32(1)
    end
    fac_ptr = Vector{Int32}(undef, n + 1)
    fac_ptr[1] = 1
    @inbounds for m in 1:n
        fac_ptr[m + 1] = fac_ptr[m] + depth[m]
    end
    fac = Vector{Int32}(undef, Int(fac_ptr[end]) - 1)
    @inbounds for m in 2:n
        pl = fac_ptr[parent[m]]
        plen = depth[parent[m]]
        base = fac_ptr[m]
        for t in 0:(plen - 1)
            fac[base + t] = fac[pl + t]
        end
        fac[base + plen] = leaf[m]
    end
    return fac, fac_ptr
end

"""Threaded monomial update: each `v[m]` is the independent product of its factor chain, so
there is no cross-iteration dependency. Bit-identical to `update_v!` (same factor order) and
to a serial run (`@batch` runs serially on one thread)."""
function update_v_flat!(v, fac, fac_ptr, u)
    Polyester.@batch for m in 2:length(v)
        s = one(ComplexF64)
        @inbounds for k in fac_ptr[m]:(fac_ptr[m + 1] - 1)
            j = fac[k]
            s *= j > 0 ? u[j] : conj(u[-j])
        end
        @inbounds v[m] = s
    end
    return v
end

@inline function _rowsum(Mt::SparseMatrixCSC, v, i)
    s = zero(eltype(v))
    @inbounds for k in Mt.colptr[i]:(Mt.colptr[i + 1] - 1)
        s += Mt.nzval[k] * v[Mt.rowval[k]]
    end
    return s
end

"""`du = M * v` via `Mt` (CSR); serial and threaded sum each row in the same order, so the
results are identical bit-for-bit regardless of `parallel`. The threaded branch uses
Polyester's persistent pool (`@batch`), which runs serially on a single thread."""
function spmv!(du, Mt::SparseMatrixCSC, v, parallel::Bool)
    if parallel
        Polyester.@batch for i in eachindex(du)
            @inbounds du[i] = _rowsum(Mt, v, i)
        end
    else
        @inbounds for i in eachindex(du)
            du[i] = _rowsum(Mt, v, i)
        end
    end
    return du
end

function (k::MomentKernel)(du, u, p, t)
    if k.parallel
        update_v_flat!(k.v, k.fac, k.fac_ptr, u)
    else
        update_v!(k.v, k.parent, k.leaf, u)
    end
    spmv!(du, k.Mt, k.v, k.parallel)
    return nothing
end

# ---- parameter sweeps without relowering ----------------------------------------------

"""Position of each `A[rows[k], cols[k]]` in `A.nzval` (duplicates accumulate into the same
slot). Called with the transposed matrix `Mt` and swapped indices, since `Mt`'s (row, col)
is (monomial, state)."""
function nz_map(A, rows, cols)
    nzmap = Vector{Int}(undef, length(rows))
    for k in eachindex(rows)
        j = cols[k]
        r = Int(A.colptr[j]):(Int(A.colptr[j + 1]) - 1)
        p = searchsortedfirst(view(A.rowval, r), rows[k]) + first(r) - 1
        @assert A.rowval[p] == rows[k]
        nzmap[k] = p
    end
    return nzmap
end

"""
The parameter payload of a kernel `ODEProblem` (`prob.p`). Carries the discovered
parameter occurrences, the current values, and everything needed to rewrite `Mt.nzval`
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

function KernelParameters(ir::MomentIR, Mt::SparseMatrixCSC, values::Dict)
    evalcoeffs = vals -> coefficient_values(ir, vals)
    return KernelParameters(
        ir.params, Dict{Any, Any}(values), evalcoeffs, ir.coo_c,
        nz_map(Mt, ir.coo_j, ir.coo_i),   # Mt is transposed: (row, col) = (monomial, state)
    )
end

function write_nzval!(k::MomentKernel, kp::KernelParameters, cvals)
    fill!(k.Mt.nzval, zero(ComplexF64))
    @inbounds for t in eachindex(kp.nzmap)
        k.Mt.nzval[kp.nzmap[t]] += cvals[kp.coo_c[t]]
    end
    return k
end
