# Compact serial numerical evaluator for MomentIR.
#
# This layer contains no symbolic metadata and no SciML integration. A MomentKernel stores
# only integer execution tables. Numerical coefficient values are supplied separately on
# each call; parameter binding belongs to the next layer.
#
# Scratch is task-local rather than kernel-local. The evaluator never yields, so two serial
# kernel calls on one task cannot overlap; kernels with the same numeric type and monomial
# count may therefore reuse one buffer. The package-global TLS key also means completed tasks
# do not retain discarded kernel objects.

const _MOMENT_SCRATCH_TLS_KEY = Ref{Nothing}()
const _MomentScratchRegistry = Dict{DataType, Any}

struct _MomentScratchCache{T}
    buffers::Dict{Int, Vector{T}}
end

_MomentScratchCache(::Type{T}) where {T} = _MomentScratchCache{T}(Dict{Int, Vector{T}}())

function _moment_scratch_registry()
    tls = task_local_storage()
    registry = get(tls, _MOMENT_SCRATCH_TLS_KEY, nothing)
    if registry === nothing
        registry = _MomentScratchRegistry()
        tls[_MOMENT_SCRATCH_TLS_KEY] = registry
    end
    return registry::_MomentScratchRegistry
end

function _moment_scratch_cache(::Type{T}) where {T}
    registry = _moment_scratch_registry()
    cache = get(registry, T, nothing)
    if cache === nothing
        cache = _MomentScratchCache(T)
        registry[T] = cache
    end
    return cache::_MomentScratchCache{T}
end

function _moment_buffer(::Type{T}, n::Int) where {T}
    cache = _moment_scratch_cache(T)
    buffer = get(cache.buffers, n, nothing)
    if buffer === nothing
        buffer = Vector{T}(undef, n)
        buffer[1] = one(T)
        cache.buffers[n] = buffer
    end
    return buffer::Vector{T}
end

"""
    MomentKernel{T}

Concrete serial execution plan compiled from a `MomentIR`. The kernel retains only integer
monomial/row tables; it does not retain symbolic states, coefficients, parameters, source
equations, or mutable scratch buffers.
"""
struct MomentKernel{T}
    parent::Vector{Int32}
    leaf::Vector{Int32}
    rowptr::Vector{Int32}
    monomial::Vector{Int32}
    coeff_id::Vector{Int32}
end

MomentKernel(ir::MomentIR) = MomentKernel(ir, ComplexF64)
function MomentKernel(ir::MomentIR, ::Type{T}) where {T <: Number}
    isconcretetype(T) || throw(ArgumentError("MomentKernel numeric type must be concrete; got $T"))
    return MomentKernel{T}(ir.parent, ir.leaf, ir.rowptr, ir.monomial, ir.coeff_id)
end

Base.length(kernel::MomentKernel) = length(kernel.rowptr) - 1

"""Evaluate every shared monomial in prefix order."""
function _update_monomials!(v::AbstractVector{T}, parent, leaf, u::AbstractVector{T}) where {T}
    @inbounds v[1] = one(T)
    @inbounds for m in 2:length(v)
        j = leaf[m]
        x = j > 0 ? u[j] : conj(u[-j])
        v[m] = v[parent[m]] * x
    end
    return v
end

"""Evaluate equation rows from shared monomials and the pooled numeric coefficient vector."""
function _evaluate_moment_rows!(
        du::AbstractVector{T},
        rowptr,
        monomial,
        coeff_id,
        coeffs::AbstractVector{T},
        v::AbstractVector{T},
    ) where {T}
    @inbounds for i in eachindex(du)
        acc = zero(T)
        for k in rowptr[i]:(rowptr[i + 1] - 1)
            acc += coeffs[coeff_id[k]] * v[monomial[k]]
        end
        du[i] = acc
    end
    return du
end

"""
    kernel(du, u, coeffs)

Evaluate the structured RHS in place. `coeffs` is the numeric value of the pooled symbolic
coefficient table in the source `MomentIR`; parameter binding and refresh are deliberately
outside this evaluator.
"""
function (kernel::MomentKernel{T})(
        du::AbstractVector{T},
        u::AbstractVector{T},
        coeffs::AbstractVector{T},
    ) where {T}
    v = _moment_buffer(T, length(kernel.parent))
    _update_monomials!(v, kernel.parent, kernel.leaf, u)
    _evaluate_moment_rows!(
        du,
        kernel.rowptr,
        kernel.monomial,
        kernel.coeff_id,
        coeffs,
        v,
    )
    return nothing
end
