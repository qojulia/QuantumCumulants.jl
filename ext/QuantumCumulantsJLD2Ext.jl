# JLD2 storage for the kernel cache (`KernelBackend(cache = path)`). The extension only
# moves bytes: the digest/verify logic and the RGF reconstruction stay in the package
# (`src/kernel_cache.jl`), so this file carries no knowledge of the IR.
module QuantumCumulantsJLD2Ext

using QuantumCumulants: QuantumCumulants
using JLD2: jldopen

function _store!(path, text, digest, payload)
    try
        jldopen(path, "a+") do f
            haskey(f, digest) && return nothing
            f["$digest/text"] = text
            for (k, v) in payload
                f["$digest/$k"] = v
            end
            return nothing
        end
    catch err
        @warn "kernel cache not writable; the kernel is built but not stored" path err
    end
    return nothing
end

"""Raw payload `Dict` or `nothing` (miss / failed byte-exact verify / unreadable file)."""
function _load(path, text, digest)
    isfile(path) || return nothing
    return try
        jldopen(path, "r") do f
            haskey(f, digest) || return nothing
            f["$digest/text"] == text || return nothing   # byte-exact verify
            Dict{String, Any}(
                k => f["$digest/$k"] for k in keys(f[digest]) if k != "text"
            )
        end
    catch err
        @warn "kernel cache unreadable; lowering fresh" path err
        nothing
    end
end

function __init__()
    QuantumCumulants._CACHE_BACKEND[] = (store! = _store!, load = _load)
    return nothing
end

end # module
