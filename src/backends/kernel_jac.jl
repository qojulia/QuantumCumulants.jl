# Analytic sparse Jacobian from the MomentIR (issue #294, M·v design).
#
# For a term c · Π_f u_{j_f} on equation i, the derivative wrt u_j is
# c · mult_j · (the monomial with one factor j removed). The delete-one complement of a
# sorted factor tuple is generally NOT a prefix, so complements are registered as
# additional monomials (bringing their own prefixes) before the kernel is materialized;
# `jacobian_ir` therefore returns an EXTENDED MomentIR to build the MomentKernel from, so
# RHS and Jacobian share one monomial id space.
#
# Holomorphic case only: a conj(u) factor has zero holomorphic derivative but a nonzero
# Wirtinger derivative ∂f/∂ū, so for conj-folded systems this J would be silently wrong.
# `jacobian_ir` throws `HolomorphicJacobianError` on any conj factor; unfolded systems
# (get_adjoints=true, the default) contain none.

struct JacIR
    Jproto::SparseMatrixCSC{ComplexF64, Int32}  # sparsity pattern, values overwritten per call
    nzptr::Vector{Int32}                        # per structural nz: range into the entry lists
    e_cid::Vector{Int32}                        # entry: pooled coefficient id
    e_mono::Vector{Int32}                       # entry: complement monomial id
    e_mult::Vector{Int32}                       # entry: multiplicity of the removed factor
end

"""Factor tuple of each monomial id, reconstructed from the prefix chains."""
function ir_factors(ir::MomentIR)
    fs = Vector{Vector{Int32}}(undef, length(ir.parent))
    fs[1] = Int32[]
    for m in 2:length(ir.parent)
        fs[m] = vcat(fs[ir.parent[m]], ir.leaf[m])   # parents precede children
    end
    return fs
end

"""
    jacobian_ir(ir) -> (ir_extended, jac_ir)

Build the Jacobian tables. `ir_extended` supersedes `ir` (same COO, possibly more
monomials); construct the `MomentKernel` from it so `v` covers the complements.
"""
function jacobian_ir(ir::MomentIR)
    factors = ir_factors(ir)
    mono_ids = Dict{Vector{Int32}, Int32}(f => Int32(m) for (m, f) in enumerate(factors))
    parent = copy(ir.parent)
    leaf = copy(ir.leaf)
    function mono_id!(fs::Vector{Int32})
        return get!(mono_ids, fs) do
            p = mono_id!(fs[1:(end - 1)])
            push!(parent, p)
            push!(leaf, fs[end])
            Int32(length(parent))
        end
    end
    # accumulate entries per Jacobian position (i, j)
    entries = Dict{Tuple{Int32, Int32}, Vector{NTuple{3, Int32}}}()
    for t in eachindex(ir.coo_i)
        i, m, cid = ir.coo_i[t], ir.coo_j[t], ir.coo_c[t]
        fs = factors[m]
        for j in unique(fs)
            j < 0 && throw(HolomorphicJacobianError())
            mult = Int32(count(==(j), fs))
            comp = fs[1:end]
            deleteat!(comp, findfirst(==(j), comp))
            push!(get!(entries, (i, j), NTuple{3, Int32}[]), (cid, mono_id!(comp), mult))
        end
    end
    ir_ext = MomentIR(
        ir.nstates, parent, leaf, ir.coeffs, ir.coo_i, ir.coo_j, ir.coo_c, ir.params
    )
    # CSC order: sort positions column-major, flatten entry lists with pointers
    pos = sort!(collect(keys(entries)); by = p -> (p[2], p[1]))
    nzptr = Int32[1]
    e_cid = Int32[]
    e_mono = Int32[]
    e_mult = Int32[]
    for p in pos
        for (cid, comp, mult) in entries[p]
            push!(e_cid, cid)
            push!(e_mono, comp)
            push!(e_mult, mult)
        end
        push!(nzptr, Int32(length(e_cid) + 1))
    end
    Jproto = sparse(
        [p[1] for p in pos], [p[2] for p in pos], zeros(ComplexF64, length(pos)),
        ir.nstates, ir.nstates,
    )
    return ir_ext, JacIR(Jproto, nzptr, e_cid, e_mono, e_mult)
end

"""Jacobian callable: fills `Jmat.nzval` (Jmat must share Jproto's sparsity pattern). `v` is
the same per-thread scratch set as the RHS kernel (shared within a `prob`, thread-local
across concurrent callers), so the Jacobian is reentrant on the same footing as the RHS."""
struct JacKernel
    jac::JacIR
    parent::Vector{Int32}
    leaf::Vector{Int32}
    c::Vector{ComplexF64}
    v::Vector{Vector{ComplexF64}}
end

function JacKernel(ir_ext::MomentIR, jac::JacIR, cvals::Vector{ComplexF64})
    return JacKernel(jac, ir_ext.parent, ir_ext.leaf, cvals, _make_vbufs(length(ir_ext.parent)))
end

function (jk::JacKernel)(Jmat, u, p, t)
    v = _vbuf(jk.v)                           # this thread's scratch (reentrant)
    update_v!(v, jk.parent, jk.leaf, u)
    nzv = Jmat.nzval
    @inbounds for k in eachindex(nzv)
        acc = zero(ComplexF64)
        for e in jk.jac.nzptr[k]:(jk.jac.nzptr[k + 1] - 1)
            acc += jk.jac.e_mult[e] * jk.c[jk.jac.e_cid[e]] * v[jk.jac.e_mono[e]]
        end
        nzv[k] = acc
    end
    return nothing
end
