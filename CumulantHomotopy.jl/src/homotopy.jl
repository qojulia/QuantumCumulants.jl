# Solve a real stationary polynomial system with HomotopyContinuation.jl (ROADMAP Phase 0,
# step 4). All complex roots are enumerated with the default (polyhedral) start system; the
# physical steady states are the real roots.

# Structural walk of a real Symbolics expression into a HomotopyContinuation `Expression`.
# Real variables map to HC variables (`varmap`); parameters are substituted with their
# numeric values (`psub`); no `Symbolics.substitute` is used.
function _to_hc(ex, varmap, psub)
    ex = _uw(ex)
    ex isa Number && return ex
    haskey(varmap, ex) && return varmap[ex]
    SymbolicUtils.isconst(ex) && return ex.val
    if !SymbolicUtils.iscall(ex)
        SymbolicUtils.issym(ex) || error("unexpected leaf $ex :: $(typeof(ex))")
        haskey(psub, ex) && return psub[ex]
        error("no numeric value supplied for parameter $ex")
    end
    op = SymbolicUtils.operation(ex)
    return op([_to_hc(a, varmap, psub) for a in SymbolicUtils.arguments(ex)]...)
end

"""
    StationaryResult

Result of enumerating stationary states of a hierarchy.

# Fields
- `moments`: one `OrderedDict(⟨O⟩ => value)` per solution, giving every hierarchy moment.
- `coordinates`: the raw real-coordinate root vectors (aligned with `system.variables`).
- `system`: the [`StationaryPolynomialSystem`](@ref) that was solved.
- `raw`: the underlying `HomotopyContinuation` result.
"""
struct StationaryResult
    moments::Vector{OrderedDict{Any,ComplexF64}}
    coordinates::Vector{Vector{ComplexF64}}
    system::StationaryPolynomialSystem
    raw::Any
end

Base.length(r::StationaryResult) = length(r.moments)
Base.getindex(r::StationaryResult, i) = r.moments[i]
Base.iterate(r::StationaryResult, s = 1) = s > length(r) ? nothing : (r.moments[s], s + 1)

function Base.show(io::IO, r::StationaryResult)
    println(io, "StationaryResult: ", length(r), " stationary state(s)")
    print(io, "Moments per state: ", length(r.system.moment_map))
end

# Reconstruct every hierarchy moment from a real-coordinate root.
function _reconstruct(system::StationaryPolynomialSystem, coord::AbstractVector)
    sub = Dict(_uw(system.variables[i]) => coord[i] for i in eachindex(system.variables))
    out = OrderedDict{Any,ComplexF64}()
    for (moment, expr) in system.moment_map
        rev = Symbolics.value(Symbolics.substitute(real(expr), sub))
        imv = Symbolics.value(Symbolics.substitute(imag(expr), sub))
        out[moment] = ComplexF64(rev) + im * ComplexF64(imv)
    end
    return out
end

"""
    stationary_states(system::StationaryPolynomialSystem, parameters; only_physical=true)
    stationary_states(eqs::MeanfieldEquations, parameters; kwargs...)

Enumerate the stationary states of a QuantumCumulants hierarchy with
HomotopyContinuation.jl. `parameters` maps each symbolic parameter to a numeric value.

With `only_physical=true` (default) only the real roots are returned; these are the physical
steady states, since the conjugation relation `⟨O†⟩ = ⟨O⟩*` is built into the real
coordinates. With `only_physical=false` every complex root is returned. Extra keyword
arguments are forwarded to `HomotopyContinuation.solve`.
"""
function stationary_states(
    system::StationaryPolynomialSystem, parameters; only_physical = true, kwargs...
)
    psub = Dict{Any,Any}(_uw(k) => v for (k, v) in pairs(parameters))
    missing_params = [p for p in system.parameters if !haskey(psub, _uw(p))]
    isempty(missing_params) ||
        error("no numeric value supplied for parameter(s): $(join(missing_params, ", "))")

    hcvars = HC.variables(:x, 1:length(system.variables))
    varmap = Dict(_uw(system.variables[i]) => hcvars[i] for i in eachindex(system.variables))
    exprs = [_to_hc(eq, varmap, psub) for eq in system.equations]
    hc_system = HC.System(exprs; variables = hcvars)

    result = HC.solve(hc_system; kwargs...)
    coords = only_physical ? HC.real_solutions(result) : HC.solutions(result)
    coords = [ComplexF64.(c) for c in coords]

    moments = [_reconstruct(system, c) for c in coords]
    return StationaryResult(moments, coords, system, result)
end

function stationary_states(eqs::MeanfieldEquations, parameters; only_physical = true, kwargs...)
    return stationary_states(realify(eqs), parameters; only_physical, kwargs...)
end
