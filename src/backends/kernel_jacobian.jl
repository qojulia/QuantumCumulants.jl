# Analytic sparse Jacobian for the direct moment backend.
#
# Jacobian construction consumes MomentIR but never mutates or extends it. Delete-one
# monomial complements live in a derivative-local prefix table, so enabling a Jacobian does
# not change the RHS kernel's monomial count or execution cost.

struct HolomorphicJacobianError <: Exception end
function Base.showerror(io::IO, ::HolomorphicJacobianError)
    return print(
        io,
        "HolomorphicJacobianError: this hierarchy contains conjugated state factors, so a " *
            "single complex n×n Jacobian is not the full derivative. Use an unfolded " *
            "closure (`get_adjoints = true`) or a solver path that does not request the " *
            "direct analytic Jacobian.",
    )
end

"""Reconstruct the canonical signed factor tuple of every MomentIR monomial."""
function _moment_ir_factor_tuples(ir::MomentIR)
    factors = Vector{Tuple}(undef, length(ir.parent))
    factors[1] = ()
    @inbounds for m in 2:length(ir.parent)
        factors[m] = (factors[ir.parent[m]]..., ir.leaf[m])
    end
    return factors
end

struct MomentJacobianKernel{T, M}
    prototype::M
    parent::Vector{Int32}
    leaf::Vector{Int32}
    nzptr::Vector{Int32}
    coeff_id::Vector{Int32}
    monomial::Vector{Int32}
    multiplicity::Vector{Int32}
end

"""
    MomentJacobianKernel(ir, T = ComplexF64)

Compile the holomorphic derivative of `ir` into a sparse structural plan. Each structural
Jacobian nonzero owns one contiguous range of derivative terms. Complement monomials are
shared within this derivative plan but remain independent of the RHS monomial table.
"""
MomentJacobianKernel(ir::MomentIR) = MomentJacobianKernel(ir, ComplexF64)
function MomentJacobianKernel(ir::MomentIR, ::Type{T}) where {T <: Number}
    isconcretetype(T) || throw(
        ArgumentError("MomentJacobianKernel numeric type must be concrete; got $T"),
    )
    any(j -> j < 0, ir.leaf) && throw(HolomorphicJacobianError())

    factors = _moment_ir_factor_tuples(ir)

    # Independent prefix table for delete-one complements.
    edges = Dict{Tuple{Int32, Int32}, Int32}()
    parent = Int32[0]
    leaf = Int32[0]
    function complement_id!(fs::Tuple)
        p = Int32(1)
        @inbounds for factor in fs
            f = Int32(factor)
            key = (p, f)
            p = get!(edges, key) do
                push!(parent, p)
                push!(leaf, f)
                Int32(length(parent))
            end
        end
        return p
    end

    # One Jacobian location can receive contributions from several RHS monomials.
    entries = Dict{Tuple{Int32, Int32}, Vector{NTuple{3, Int32}}}()
    for i in 1:length(ir)
        for term in ir.rowptr[i]:(ir.rowptr[i + 1] - 1)
            fs = factors[ir.monomial[term]]
            q = 1
            while q <= length(fs)
                j = Int32(fs[q])
                j < 0 && throw(HolomorphicJacobianError())
                r = q + 1
                while r <= length(fs) && fs[r] == j
                    r += 1
                end
                mult = Int32(r - q)
                complement = (fs[1:(q - 1)]..., fs[(q + 1):end]...)
                mid = complement_id!(complement)
                push!(
                    get!(entries, (Int32(i), j), NTuple{3, Int32}[]),
                    (ir.coeff_id[term], mid, mult),
                )
                q = r
            end
        end
    end

    # Sparse CSC order is column-major; flatten contribution lists in that same order.
    positions = sort!(collect(keys(entries)); by = p -> (p[2], p[1]))
    rows = Int32[p[1] for p in positions]
    cols = Int32[p[2] for p in positions]
    prototype = sparse(rows, cols, zeros(T, length(positions)), length(ir), length(ir))

    nzptr = Int32[1]
    coeff_id = Int32[]
    monomial = Int32[]
    multiplicity = Int32[]
    for position in positions
        for (cid, mid, mult) in entries[position]
            push!(coeff_id, cid)
            push!(monomial, mid)
            push!(multiplicity, mult)
        end
        push!(nzptr, Int32(length(coeff_id) + 1))
    end

    return MomentJacobianKernel{T, typeof(prototype)}(
        prototype,
        parent,
        leaf,
        nzptr,
        coeff_id,
        monomial,
        multiplicity,
    )
end

function (kernel::MomentJacobianKernel{T})(
        J::SparseMatrixCSC{T},
        u::AbstractVector{T},
        p::KernelParameters{T},
        t,
    ) where {T}
    v = _moment_buffer(T, length(kernel.parent))
    _update_monomials!(v, kernel.parent, kernel.leaf, u)
    nzval = J.nzval
    @inbounds for k in eachindex(nzval)
        acc = zero(T)
        for entry in kernel.nzptr[k]:(kernel.nzptr[k + 1] - 1)
            acc +=
                kernel.multiplicity[entry] *
                p.coeffs[kernel.coeff_id[entry]] *
                v[kernel.monomial[entry]]
        end
        nzval[k] = acc
    end
    return nothing
end
