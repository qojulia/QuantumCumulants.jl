# Power spectrum of a `CorrelationFunction`: the Laplace-transformed two-time correlation
# and the compiled numeric kernel that solves the τ-system linear response. Included after
# correlation.jl, so it builds on the shared steady-state resolver and ambient helpers
# (`_ss_resolver`, `_undo_ancilla`, `_ambient_param`, `_avg_conj_of`, `_scalarize`).

# Compiled numeric kernel for the spectrum linear system, built once per `Spectrum`
# from the symbolic τ-RHSs (`_build_spectrum_kernel`) and cached. `f_Alin`/`f_Aalin`/`f_b`
# are `build_function`-compiled closures taking the per-call input vector (params, τ,
# ambient leaf values) in the order implied by `param_syms`, then τ, `ambient_leaves`,
# `extra_syms`. `f_Aalin` is `nothing` when the system has no anti-linear coupling.
struct SpectrumKernel{FL, FA, FB}
    n::Int
    has_alin::Bool
    f_Alin::FL
    f_Aalin::FA
    f_b::FB
    param_syms::Vector{Any}
    ambient_leaves::Vector{Any}
    extra_syms::Vector{Any}
end

"""
    Spectrum(c::CorrelationFunction, ps)

The Laplace-transformed correlation as a callable: `S(ω, u_end, p0)` returns the
symmetric power spectrum `2·Re{∫₀^∞ g(τ)e^{-iωτ}dτ}`. Solves `(iω·I - A)X̃ = x̃(0)`
where `A` is the τ-system Jacobian at steady state and `x̃ = x - x_∞` is centred so
the one-sided transform converges. The linear/anti-linear split (the conjugate
columns from `get_adjoints=false`) is handled by the 2n augmentation.

# Examples
```jldoctest
julia> h = FockSpace(:cavity);

julia> @qnumbers a::Destroy(h);

julia> eqs = meanfield([a' * a], a' * a, [a]; rates = [1.0], order = 2);

julia> c = CorrelationFunction(a', a, eqs);

julia> Spectrum(c)
ℱ(⟨a' * a⟩)(ω)
```
"""
struct Spectrum{C <: CorrelationFunction, P}
    c::C
    ω::Symbolics.Num
    ps::P
    # Lazily-built compiled kernel; depends only on `c.eqs` (fixed for the Spectrum's
    # lifetime), so it is built once on first evaluation and reused. Ignored by equality
    # and printing.
    cache::Base.RefValue{Union{Nothing, SpectrumKernel}}
end
function Spectrum(c::CorrelationFunction, ps)
    c.steady_state || throw(
        ArgumentError(
            "Spectrum requires `CorrelationFunction(...; steady_state = true)`. " *
                "Got `steady_state = false`."
        ),
    )
    return Spectrum{typeof(c), typeof(ps)}(
        c, first(@variables ω_spectrum), ps,
        Base.RefValue{Union{Nothing, SpectrumKernel}}(nothing),
    )
end
Spectrum(c::CorrelationFunction) = Spectrum(c, ())

function (S::Spectrum)(ω_vals::AbstractVector, u_end, p0)
    A, rhs_b, n = _spectrum_kernel(S, u_end, p0)
    I_n = Matrix{ComplexF64}(I, n, n)
    out = Vector{Float64}(undef, length(ω_vals))
    @inbounds for (k, ω_val) in enumerate(ω_vals)
        x = ((im * ω_val) .* I_n .- A) \ rhs_b
        out[k] = 2 * real(x[1])
    end
    return out
end
(S::Spectrum)(ω_val::Real, u_end, p0) = first(S([ω_val], u_end, p0))

"""
Assemble the linear system `(iω·I - A)X̃ = x̃(0)` for the spectrum. The τ-system is
linear, so its Jacobian `A` and constant term `b` are obtained from a compiled,
cached kernel (`_build_spectrum_kernel`): each RHS leaf is classified once as linear
(matches a state), anti-linear (conjugate of a state), or ambient (a steady-state
coefficient); the linear/anti-linear coefficient blocks are differentiated symbolically
and `build_function`-compiled, and the 2n augmentation closes the conjugate response.
The build runs once per `Spectrum`; this call only evaluates numerically and solves.
"""
function _spectrum_kernel(S::Spectrum, u_end, p0)
    c = S.c
    eqs = c.eqs
    kern = S.cache[]
    if kern === nothing
        kern = _build_spectrum_kernel(S)
        S.cache[] = kern
    end
    n = kern.n
    resolve = _ss_resolver(c, u_end)
    A_lin, A_alin, b_const = _eval_spectrum_blocks(kern, S, u_end, p0, resolve)
    u_τ = ComplexF64[
        resolve(_undo_ancilla(c, SymbolicUtils.unwrap(s))) for s in eqs.states
    ]
    if !kern.has_alin
        rhs_b = any(!iszero, b_const) ? u_τ .+ (A_lin \ b_const) : u_τ
        return A_lin, rhs_b, n
    end
    M = [A_lin A_alin; conj.(A_alin) conj.(A_lin)]
    b_aug = vcat(b_const, conj.(b_const))
    u_aug = vcat(u_τ, conj.(u_τ))
    rhs_b = any(!iszero, b_aug) ? u_aug .+ (M \ b_aug) : u_aug
    return M, rhs_b, 2n
end

# Collect the distinct unwrapped free symbols across all entries of a symbolic block.
function _block_free_syms(block)
    out = SymbolicUtils.BasicSymbolic[]
    seen = Set{SymbolicUtils.BasicSymbolic}()
    for e in block
        eu = SymbolicUtils.unwrap(e)
        eu isa SymbolicUtils.BasicSymbolic || continue
        for v in Symbolics.get_variables(eu)
            vu = SymbolicUtils.unwrap(v)
            vu isa SymbolicUtils.BasicSymbolic || continue
            vu in seen || (push!(seen, vu); push!(out, vu))
        end
    end
    return out
end

# Error if a state placeholder survives into a coefficient block: that means a term is
# bilinear in two state moments, i.e. the τ-system is not linear (a violated contract the
# old finite-difference probe would have silently dropped).
function _assert_linear_blocks(blocks, state_ph_set)
    for block in blocks
        for v in _block_free_syms(block)
            if v in state_ph_set
                throw(
                    ArgumentError(
                        "spectrum τ-system is not linear in its states (a coefficient depends " *
                            "on a state moment). The quantum regression theorem requires a linear " *
                            "system; check the cumulant order and closure."
                    )
                )
            end
        end
    end
    return nothing
end

"""
Build the cached numeric spectrum kernel for `S`. Placeholder-resolves every RHS leaf to a
plain `Number` symbol (so `Symbolics.jacobian`/`build_function` do not descend into the
`average` envelope), classifies leaves into linear/anti-linear state directions and ambient
coefficients, differentiates the coefficient blocks `A_lin`/`A_alin` symbolically, takes the
constant term `b` by zeroing the state placeholders, and `build_function`-compiles the three
blocks over an explicit input set (`S.ps` params, `τ`, ambient leaves). `cse` toggles
common-subexpression elimination (off by default; on works but is unbenchmarked).
"""
function _build_spectrum_kernel(S::Spectrum; cse::Bool = false)
    c = S.c
    eqs = c.eqs
    n = length(eqs)
    rhss_u = [SymbolicUtils.unwrap(eq.rhs) for eq in eqs.equations]

    # Distinct RHS leaves (multiset-deduped by identity), same as the old kernel.
    rhs_leaves = SymbolicUtils.BasicSymbolic[]
    seen = Set{SymbolicUtils.BasicSymbolic}()
    for rhs in rhss_u, l in eachleaf(rhs)
        l in seen || (push!(seen, l); push!(rhs_leaves, l))
    end

    # Placeholder-resolve ALL leaves to plain Number symbols. jacobian/cse crash on raw
    # average leaves (`MethodError(avg, ...)`); plain symbols sidestep that.
    phs = Any[
        SymbolicUtils.unwrap(Symbolics.variable(Symbol("__spec_x", i); T = Number))
            for i in eachindex(rhs_leaves)
    ]
    leaf_to_ph = Dict{Any, Any}(rhs_leaves[i] => phs[i] for i in eachindex(rhs_leaves))
    rhss_plain = [SymbolicUtils.substitute(rhs, leaf_to_ph) for rhs in rhss_u]

    # Classify each leaf by Hermitian conjugation (same side = linear, opposite = anti-linear)
    # via the conjugate representative, robust for scaled states where adjoint does not
    # round-trip. Unmatched leaves are ambient coefficients (numeric per call).
    spec_ctx = build_ctx(eqs)
    state_map = MomentMap(
        spec_ctx, _treatments(eqs, spec_ctx),
        [undo_average(s) for s in eqs.states],
        collect(1:n),
    )
    lin_cols = Dict{Int, Vector{Int}}()
    alin_cols = Dict{Int, Vector{Int}}()
    classified = Set{Int}()
    for (k, leaf) in enumerate(rhs_leaves)
        r = match_moment(state_map, undo_average(leaf))
        r === nothing && continue
        i, same = r
        push!(get!(same ? lin_cols : alin_cols, i, Int[]), k)
        push!(classified, k)
    end
    ambient_cols = [k for k in eachindex(rhs_leaves) if !(k in classified)]

    # One full-vector jacobian over all state placeholders; per-state subset calls are
    # avoided (only the full-vector form reproduces the cross terms correctly).
    statecols = sort(unique(vcat(values(lin_cols)..., values(alin_cols)..., Int[])))
    statephs = Any[phs[k] for k in statecols]
    colpos = Dict(statecols[i] => i for i in eachindex(statecols))
    Jfull = isempty(statephs) ? Matrix{Any}(undef, n, 0) :
        Symbolics.jacobian(rhss_plain, statephs)
    # A_block[:, j] = sum of jacobian columns of state j's leaves (matches the old probe,
    # which set all of a state's leaves to 1 simultaneously).
    assemble = function (cols)
        A = Matrix{Any}(undef, n, n)
        fill!(A, 0)
        for (j, ks) in cols, kk in ks, i in 1:n
            A[i, j] = SymbolicUtils.unwrap(A[i, j]) + SymbolicUtils.unwrap(Jfull[i, colpos[kk]])
        end
        return A
    end
    Alin = assemble(lin_cols)
    has_alin = !isempty(alin_cols)
    Aalin = has_alin ? assemble(alin_cols) : nothing
    zero_states = Dict{Any, Any}(phs[k] => 0 for k in statecols)
    bvec = [SymbolicUtils.substitute(rhs, zero_states) for rhs in rhss_plain]

    state_ph_set = Set(statephs)
    blocks = has_alin ? Any[Alin, Aalin, bvec] : Any[Alin, bvec]
    _assert_linear_blocks(blocks, state_ph_set)

    # Explicit input set: params (S.ps), τ, ambient placeholders, plus any stray params
    # appearing in the blocks but not passed to `Spectrum` (matched to 0 per call, as the
    # old kernel did via `_scalarize`).
    param_syms = Any[SymbolicUtils.unwrap(p) for p in S.ps]
    ambient_phs = Any[phs[k] for k in ambient_cols]
    known = Set{SymbolicUtils.BasicSymbolic}()
    for s in param_syms
        push!(known, SymbolicUtils.unwrap(s))
    end
    for s in statephs
        push!(known, SymbolicUtils.unwrap(s))
    end
    for s in ambient_phs
        push!(known, SymbolicUtils.unwrap(s))
    end
    push!(known, SymbolicUtils.unwrap(c.τ))
    # `get_variables` reports the imaginary unit `IM` as a free symbol; it is a constant
    # emitted literally by `build_function`, never an input (binding it would zero `x*im`).
    push!(known, SymbolicUtils.unwrap(Symbolics.IM))
    extra_syms = Any[]
    for block in blocks, v in _block_free_syms(block)
        (v in known) && continue
        push!(known, v)
        push!(extra_syms, v)
    end
    ambient_leaves = Any[rhs_leaves[k] for k in ambient_cols]
    inputs = vcat(
        param_syms, Any[SymbolicUtils.unwrap(c.τ)], ambient_phs, extra_syms,
    )
    compile = function (M)
        ex = Symbolics.build_function(M, inputs; expression = Val{false}, cse = cse)
        return ex isa Tuple ? ex[1] : ex
    end
    f_Alin = compile(Alin)
    f_Aalin = has_alin ? compile(Aalin) : nothing
    f_b = compile(bvec)
    return SpectrumKernel(
        n, has_alin, f_Alin, f_Aalin, f_b, param_syms, ambient_leaves, extra_syms,
    )
end

# Evaluate the cached kernel's coefficient blocks numerically for `(u_end, p0)`. Builds the
# parameter/ambient substitution exactly as the old kernel, then calls the compiled blocks
# with the input vector in the order fixed at build time.
function _eval_spectrum_blocks(kern::SpectrumKernel, S::Spectrum, u_end, p0, resolve)
    c = S.c
    eqs = c.eqs
    p_sub = _build_p_sub(S.ps, p0, c.τ, u_end, eqs)
    for avg in c.ambient
        aons = SQA.acts_on(avg)
        lookup = (c.aon_ancilla in aons && length(aons) == 1) ? _undo_ancilla(c, avg) : avg
        p_sub[avg] = resolve(lookup)
    end
    invals = ComplexF64[]
    for p in kern.param_syms
        push!(invals, ComplexF64(get(p_sub, p, 0)))
    end
    push!(invals, 0.0 + 0.0im)  # τ
    for leaf in kern.ambient_leaves
        push!(invals, ComplexF64(get(p_sub, leaf, 0)))
    end
    for s in kern.extra_syms
        push!(invals, ComplexF64(get(p_sub, s, 0)))
    end
    A_lin = ComplexF64.(kern.f_Alin(invals))
    A_alin = kern.has_alin ? ComplexF64.(kern.f_Aalin(invals)) : nothing
    b_const = ComplexF64.(kern.f_b(invals))
    return A_lin, A_alin, b_const
end

function _build_p_sub(ps, p0, τ, u_end, eqs)
    sub = Dict{Any, Any}()
    if ps !== nothing && !isempty(ps)
        for (p, v) in zip(ps, p0)
            sub[SymbolicUtils.unwrap(p)] = ComplexF64(v)
        end
    end
    if u_end isa AbstractDict
        for (k, v) in u_end
            ku = SymbolicUtils.unwrap(k)
            sub[ku] = ComplexF64(v)
            SQA.is_average(ku) && (sub[SymbolicUtils.unwrap(_avg_conj_of(ku))] = ComplexF64(conj(v)))
        end
    elseif u_end isa AbstractVector
        for (avg, v) in zip(eqs.states, u_end)
            sub[SymbolicUtils.unwrap(avg)] = ComplexF64(v)
        end
    end
    sub[SymbolicUtils.unwrap(τ)] = 0.0
    return sub
end
