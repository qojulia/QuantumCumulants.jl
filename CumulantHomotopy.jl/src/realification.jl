# Realification: turn a completed QuantumCumulants stationary hierarchy into a minimal
# *real* polynomial system, as specified in SPEC §3.2.
#
# Each complex moment ⟨O⟩ is written in terms of real coordinates,
#
#     ⟨O⟩  = Oᵣ + i·Oᵢ ,   ⟨O†⟩ = Oᵣ - i·Oᵢ      (non-Hermitian conjugate pair)
#     ⟨O⟩  = Oᵣ                                    (Hermitian moment)
#
# The stationary equation of one representative per conjugate class is split into its real
# and imaginary parts. The redundant conjugate equation is dropped. The result is a square
# real polynomial system whose *real* roots are the physical steady states (the conjugation
# relation ⟨O†⟩ = ⟨O⟩* is built into the coordinates rather than checked afterwards).

const _uw = SymbolicUtils.unwrap

# Is `v` an ⟨…⟩ average leaf?
_is_average(v) = SymbolicUtils.iscall(v) && SymbolicUtils.operation(v) isa SQA.AvgFunc

# The conjugate moment leaf ⟨O†⟩ of an average leaf ⟨O⟩.
_conj_leaf(v) = _uw(SQA.average(SQA.undo_average(v)'))

# A concrete zero test for an expanded symbolic expression (Symbolics.iszero is symbolic).
_iszero(x) = (v = _uw(x); v isa Number ? iszero(v) : (SymbolicUtils.isconst(v) && iszero(v.val)))

# A valid-identifier stem for the real variables of a moment, e.g. ⟨a' * a⟩ → "a⁺a".
function _moment_stem(v)
    return replace(string(v), "⟨" => "", "⟩" => "", "*" => "", " " => "", "'" => "⁺")
end

"""
    StationaryPolynomialSystem

A minimal real polynomial system exported from a completed QuantumCumulants hierarchy
(SPEC §3.2). Its real roots are the stationary states.

# Fields
- `equations`: real polynomial residuals, each understood as `== 0`.
- `variables`: the real coordinates `Oᵣ`, `Oᵢ`.
- `parameters`: the remaining symbolic parameters (to be given numeric values when solving).
- `moment_map`: ordered map from each hierarchy moment `⟨O⟩` to its real-coordinate
  expression `Oᵣ + i·Oᵢ` (a `Complex{Num}`), used to reconstruct complex moments from a root.
"""
struct StationaryPolynomialSystem
    equations::Vector{Symbolics.Num}
    variables::Vector{Symbolics.Num}
    parameters::Vector{Symbolics.Num}
    moment_map::OrderedDict{Any,Complex{Symbolics.Num}}
end

function Base.show(io::IO, s::StationaryPolynomialSystem)
    println(io, length(s.equations), " real stationary equations")
    println(io, "Variables: ", join(string.(s.variables), ", "))
    print(io, "Parameters: ", join(string.(s.parameters), ", "))
end

# Rebuild `ex` replacing each average leaf by its real-coordinate expression (`leafsub`),
# QuantumCumulants' imaginary unit `Symbolics.IM` by the genuine complex `im` (so that
# `real`/`imag` fold `im^2 → -1`), and leaving parameters symbolic. A structural walk,
# never `Symbolics.substitute`, which trips on the average leaves' shape metadata.
function _rebuild_real(ex, leafsub)
    ex = _uw(ex)
    ex isa Number && return ex
    haskey(leafsub, ex) && return leafsub[ex]
    ex === Symbolics.IM && return im
    SymbolicUtils.isconst(ex) && return ex.val
    SymbolicUtils.iscall(ex) || return Symbolics.Num(ex)
    op = SymbolicUtils.operation(ex)
    return op([_rebuild_real(a, leafsub) for a in SymbolicUtils.arguments(ex)]...)
end

"""
    realify(eqs::MeanfieldEquations) -> StationaryPolynomialSystem

Convert a completed (closed) QuantumCumulants stationary hierarchy into a real polynomial
system (SPEC §3.2). `eqs` must be closed: every moment on a right-hand side is a state or
the conjugate of one. Call [`complete`](@ref) first otherwise.
"""
function realify(eqs::MeanfieldEquations)
    states = _uw.(eqs.states)
    state_eq = Dict(states[i] => _uw(eqs.equations[i].rhs) for i in eachindex(states))

    # Conjugation-closed set of moment leaves (states plus any conjugate leaves on the RHS).
    leaves = Set{Any}()
    for rhs in values(state_eq), v in Symbolics.get_variables(rhs)
        _is_average(v) && push!(leaves, v)
    end
    for s in states
        push!(leaves, s)
    end
    for v in collect(leaves)
        push!(leaves, _conj_leaf(v))
    end
    ordered_leaves = sort(collect(leaves); by = string)

    # Assign real coordinates per conjugate class and record the representative leaves.
    leafsub = Dict{Any,Any}()
    variables = Symbolics.Num[]
    representatives = Any[]
    used_names = Set{Symbol}()
    seen = Set{Any}()
    declare(stem, suffix) = begin
        name = Symbol(stem, suffix)
        while name in used_names
            name = Symbol(stem, "_", length(used_names), suffix)
        end
        push!(used_names, name)
        (Symbolics.@variables $name::Real)[1]
    end
    for v in ordered_leaves
        v in seen && continue
        cv = _conj_leaf(v)
        if isequal(v, cv)                                   # Hermitian moment → one real var
            vr = declare(_moment_stem(v), "ᵣ")
            push!(variables, vr)
            leafsub[v] = vr
            push!(representatives, v)
            push!(seen, v)
        else                                                # conjugate pair → (Oᵣ, Oᵢ)
            rep = haskey(state_eq, v) ? v : (haskey(state_eq, cv) ? cv : v)
            other = isequal(rep, v) ? cv : v
            stem = _moment_stem(rep)
            vr = declare(stem, "ᵣ")
            vi = declare(stem, "ᵢ")
            push!(variables, vr)
            push!(variables, vi)
            leafsub[rep] = vr + im * vi
            leafsub[other] = vr - im * vi
            push!(representatives, rep)
            push!(seen, rep)
            push!(seen, other)
        end
    end

    # One complex equation per class → real and imaginary parts.
    equations = Symbolics.Num[]
    for rep in representatives
        rhs = haskey(state_eq, rep) ? state_eq[rep] : SQA.inner_adjoint(state_eq[_conj_leaf(rep)])
        z = Symbolics.expand(_rebuild_real(rhs, leafsub))
        re = Symbolics.expand(Symbolics.real(z))
        impart = Symbolics.expand(Symbolics.imag(z))
        _iszero(re) || push!(equations, re)
        _iszero(impart) || push!(equations, impart)
    end

    if length(equations) != length(variables)
        error(
            "realification produced $(length(equations)) equations for $(length(variables)) " *
            "variables. The hierarchy may be unclosed; call `complete(eqs)` first.",
        )
    end

    # Parameters: symbolic leaves of the equations that are not real variables.
    varset = Set(_uw.(variables))
    params = Set{Any}()
    for eq in equations, v in Symbolics.get_variables(eq)
        vu = _uw(v)
        (vu in varset || vu === Symbolics.IM) && continue
        SymbolicUtils.issym(vu) && push!(params, vu)
    end
    parameters = sort(Symbolics.Num.(collect(params)); by = string)

    moment_map = OrderedDict{Any,Complex{Symbolics.Num}}()
    for s in eqs.states
        moment_map[s] = Complex{Symbolics.Num}(leafsub[_uw(s)])
    end

    return StationaryPolynomialSystem(equations, variables, parameters, moment_map)
end
