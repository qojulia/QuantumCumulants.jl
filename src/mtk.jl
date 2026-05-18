function _stable_avg_name(avg::SymbolicUtils.BasicSymbolic)
    @assert SecondQuantizedAlgebra.is_average(avg) "expected an Average BasicSymbolic"
    op = SecondQuantizedAlgebra.undo_average(avg)
    s = "avg_" * _op_name_chunk(op)
    return Symbol(s)
end

function _op_name_chunk(op::QAdd)
    isempty(op.arguments) && return "zero"
    chunks = String[]
    for (term, _) in op.arguments
        push!(chunks, join((_op_name_chunk(o) for o in term.ops), "_"))
    end
    return join(chunks, "_plus_")
end
function _op_name_chunk(op::SecondQuantizedAlgebra.QSym)
    base = string(op.name)
    extra = _op_name_extra(op)
    return isempty(extra) ? base : base * extra
end

function _op_name_extra(op::SecondQuantizedAlgebra.QSym)
    return _op_index_suffix(op)
end
function _op_name_extra(op::SecondQuantizedAlgebra.Transition)
    i = op.i isa Symbol ? string(op.i) : string(Int(op.i))
    j = op.j isa Symbol ? string(op.j) : string(Int(op.j))
    return "_" * i * j * _op_index_suffix(op)
end
_op_name_extra(op::SecondQuantizedAlgebra.Destroy) = _op_index_suffix(op)
_op_name_extra(op::SecondQuantizedAlgebra.Create) = "_dag" * _op_index_suffix(op)
function _op_name_extra(op::SecondQuantizedAlgebra.Pauli)
    return "_" * string(Int(op.axis)) * _op_index_suffix(op)
end
function _op_name_extra(op::SecondQuantizedAlgebra.Spin)
    return "_" * string(Int(op.axis)) * _op_index_suffix(op)
end

function _op_index_suffix(op::SecondQuantizedAlgebra.QSym)
    isdefined(op, :index) || return ""
    idx = op.index
    idx === SecondQuantizedAlgebra.NO_INDEX && return ""
    return "_" * string(idx.name)
end

function _avg_to_var_dict(eqs::AbstractMeanFieldEquations)
    iv   = eqs.iv
    dict = IdDict{SymbolicUtils.BasicSymbolic, Symbolics.Num}()
    dvs  = Symbolics.Num[]
    for avg in eqs.states
        v = _make_time_dependent_var(_stable_avg_name(avg), iv)
        dict[avg] = v
        push!(dvs, v)
    end
    return dict, dvs
end

function _make_time_dependent_var(name::Symbol, iv::Symbolics.Num)
    v = first(@variables $name(iv))
    return v
end

function _collect_params!(set, x, dict, iv_uw)
    if x isa SymbolicUtils.BasicSymbolic
        SymbolicUtils.isconst(x) && return
        if !SymbolicUtils.iscall(x)
            if !haskey(dict, x) && x !== iv_uw && SymbolicUtils.symtype(x) <: Real
                push!(set, x)
            end
            return
        end
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if length(args) == 1 && args[1] === iv_uw && op isa SymbolicUtils.BasicSymbolic
            return
        end
        for a in args
            _collect_params!(set, a, dict, iv_uw)
        end
    end
    return
end

"""
    to_system(eqs::MeanFieldEquations; name::Symbol)

Build a `ModelingToolkitBase.System` from the QC equation set. Substitutes
Averages with real-typed `u(t)` Num variables and passes `dvs`/`ps` explicitly.
"""
function to_system(eqs::NoiseMeanFieldEquations{O,H,Op,Jt,Jdt,R,E,S,Forward};
                   name::Symbol) where {O,H,Op,Jt,Jdt,R,E,S}
    return _to_system_sde(eqs, name, +1)
end

function to_system(eqs::NoiseMeanFieldEquations{O,H,Op,Jt,Jdt,R,E,S,Backward};
                   name::Symbol) where {O,H,Op,Jt,Jdt,R,E,S}
    return _to_system_sde(eqs, name, -1)
end

function _to_system_sde(eqs::NoiseMeanFieldEquations, name::Symbol, sign::Int)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D  = Symbolics.Differential(iv)
    dict, dvs = _avg_to_var_dict(eqs)
    drift = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs      = SymbolicUtils.substitute(eq.rhs, dict)
        drift[i] = D(dict[eq.lhs]) ~ sign * rhs
        _collect_params!(ps_set, rhs, dict, iv_uw)
    end
    ps = [MTK.toparam(p) for p in ps_set]
    return MTK.System(drift, iv, dvs, ps; name=name)
end

function to_system(eqs::MeanFieldEquations; name::Symbol)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D  = Symbolics.Differential(iv)
    dict, dvs = _avg_to_var_dict(eqs)
    conj_dict = _conj_substitution_dict(eqs, dict)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set  = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs_conj = _substitute_conj_avgs(eq.rhs, conj_dict)
        rhs      = SymbolicUtils.substitute(rhs_conj, dict)
        new_eqs[i] = D(dict[eq.lhs]) ~ rhs
        _collect_params!(ps_set, rhs, dict, iv_uw)
    end
    ps_old = collect(ps_set)
    ps = [MTK.toparam(p) for p in ps_old]
    return MTK.System(new_eqs, iv, dvs, ps; name=name)
end

# Build a substitution `⟨op†⟩ → conj(state_var(⟨op⟩))` for every state
# whose conjugate is *not* itself a state. This lets `to_system` codegen
# substitute the conjugate of a state without needing the conjugate to
# also be added to `eqs.states`. (See completion.jl::find_missing: by
# default, conjugates of states are excluded from missing-state scans.)
function _conj_substitution_dict(eqs::AbstractMeanFieldEquations,
                                 var_dict::AbstractDict)
    states = Set(eqs.states)
    conj_dict = Dict{SymbolicUtils.BasicSymbolic, Any}()
    for s in eqs.states
        cs = _avg_conj_for_codegen(s)
        cs === s && continue
        cs in states && continue
        haskey(var_dict, s) || continue
        conj_dict[cs] = conj(var_dict[s])
    end
    return conj_dict
end

function _avg_conj_for_codegen(x::SymbolicUtils.BasicSymbolic)
    SecondQuantizedAlgebra.is_average(x) || return x
    SymbolicUtils.iscall(x) || return x
    SymbolicUtils.operation(x) === SecondQuantizedAlgebra.sym_average || return x
    op = SecondQuantizedAlgebra.undo_average(x)
    return average(adjoint(op))
end

function _substitute_conj_avgs(x, conj_dict)
    isempty(conj_dict) && return x
    x isa SymbolicUtils.BasicSymbolic || return x
    if haskey(conj_dict, x)
        return conj_dict[x]
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    op === SecondQuantizedAlgebra.sym_average && return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_substitute_conj_avgs(a, conj_dict) for a in args]
    return op(new_args...)
end

"""
    initial_values(eqs::AbstractMeanFieldEquations; defaults=Dict())

Return `Dict{Symbolics.Num, ComplexF64}` mapping each state's u(t) variable to
its initial value. Unspecified averages default to `zero(ComplexF64)`.
"""
function initial_values(eqs::AbstractMeanFieldEquations;
                        defaults::AbstractDict = Dict())
    dict, _ = _avg_to_var_dict(eqs)
    u0 = Dict{Symbolics.Num, ComplexF64}()
    for avg in eqs.states
        u0[dict[avg]] = ComplexF64(get(defaults, avg, 0))
    end
    return u0
end

"""
    get_solution(sol, avg_or_op, eqs)

Query an ODESolution `sol` for the trajectory of `avg_or_op`. Accepts either a
raw `Average` BasicSymbolic or a `QField` (which is averaged internally).
"""
function get_solution(sol, avg::SymbolicUtils.BasicSymbolic,
                      eqs::AbstractMeanFieldEquations)
    dict, _ = _avg_to_var_dict(eqs)
    var = dict[avg]
    return τ -> sol(τ; idxs=var)
end
get_solution(sol, op::QField, eqs::AbstractMeanFieldEquations) =
    get_solution(sol, average(op), eqs)
