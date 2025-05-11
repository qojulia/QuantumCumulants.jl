"""
    struct CorrelationFunction

Type representing the two-time first-order correlation function of two operators.
"""
struct CorrelationFunction{OP1,OP2,OP0,DE0,DE,S}
    op1::OP1
    op2::OP2
    op2_0::OP0
    de0::DE0
    de::DE
    steady_state::S
end

"""
    CorrelationFunction(op1,op2,de0;steady_state=false,add_subscript=0,mix_choice=maximum)

The first-order two-time correlation function of two operators.

The first-order two-time correlation function of `op1` and `op2` evolving under
the system `de0`. The keyword `steady_state` determines whether the original
system `de0` was evolved up to steady state. The arguments `add_subscript`
defines the subscript added to the name of `op2` representing the constant time.

Note that the correlation function is stored in the first index of the underlying
system of equations.
"""
function CorrelationFunction(op1,op2,de0::AbstractMeanfieldEquations;
                            steady_state=false, add_subscript=0,
                            filter_func=nothing, mix_choice=maximum,
                            iv=nothing,
                            order=nothing,
                            simplify=true, kwargs...)
    if iv === nothing
        iv_ = SymbolicUtils.Sym{Real}(:τ)
        iv_ = SymbolicUtils.setmetadata(iv_, MTK.VariableSource, (:parameters, :τ))
        iv_ = SymbolicUtils.setmetadata(iv_, MTK.MTKVariableTypeCtx,MTK.PARAMETER)
    else
        iv_ = iv
    end
    h1 = hilbert(op1)
    h2 = _new_hilbert(hilbert(op2), acts_on(op2))
    h = h1⊗h2

    H0 = de0.hamiltonian
    J0 = de0.jumps
    Jd0 = de0.jumps_dagger

    op1_ = _new_operator(op1, h)
    op2_ = _new_operator(op2, h, length(h.spaces); add_subscript=add_subscript)
    op2_0 = _new_operator(op2, h)
    H = _new_operator(H0, h)
    J = [_new_operator(j, h) for j in J0]
    Jd = [_new_operator(j, h) for j in Jd0]
    lhs_new = [_new_operator(l, h) for l in de0.states]

    order_ = if order===nothing
        if de0.order===nothing
            de0.order
            order_lhs = maximum(get_order(l) for l in de0.states)
            order_corr = get_order(op1_*op2_)
            max(order_lhs, order_corr)
        else
            de0.order
        end
    else
        order
    end
    op_ = op1_*op2_
    @assert get_order(op_) <= order_

    de = meanfield(op_,H,J;Jdagger=Jd,rates=de0.rates,iv=iv_,order=order_)
    _complete_corr!(de, length(h.spaces), lhs_new, order_, steady_state;
                            filter_func=filter_func,
                            mix_choice=mix_choice,
                            simplify=simplify,
                            kwargs...)

    varmap = make_varmap(lhs_new, de0.iv)
    de0_ = begin
        eqs = Symbolics.Equation[]
        eqs_op = Symbolics.Equation[]
        ops = map(undo_average, lhs_new)
        for i=1:length(de0.equations)
            rhs = _new_operator(de0.equations[i].rhs, h)
            rhs_op = _new_operator(de0.operator_equations[i].rhs, h)
            push!(eqs, Symbolics.Equation(lhs_new[i], rhs))
            push!(eqs_op, Symbolics.Equation(ops[i], rhs_op))
        end
        MeanfieldEquations(eqs,eqs_op,lhs_new,ops,H,J,Jd,de0.rates,de0.iv,varmap,order_)
    end

    return CorrelationFunction(op1_, op2_, op2_0, de0_, de, steady_state)
end

"""
    correlation_u0(c::CorrelationFunction, u_end)

Find the vector containing the correct initial values when numerical solving
the time evolution for the correlation function.

See also: [`CorrelationFunction`](@ref) [`correlation_p0`](@ref)
"""
function correlation_u0(c::CorrelationFunction, u_end)
    a0 = c.op2_0
    a1 = c.op2
    subs = Dict(a1=>a0)
    ops = c.de.operators
    lhs = [inorder!(average(substitute(op, subs))) for op in ops]
    u0 = complex(eltype(u_end))[]
    lhs0 = c.de0.states
    τ = MTK.get_iv(c.de)
    keys = []
    for j=1:length(lhs)
        l=lhs[j]
        l_adj = _inconj(l)
        if l ∈ Set(lhs0)
            i = findfirst(isequal(l), lhs0)
            push!(u0, u_end[i])
            push!(keys, make_var(c.de.equations[j].lhs, τ))
        elseif l_adj ∈ Set(lhs0)
            i = findfirst(isequal(l_adj), lhs0)
            push!(u0, conj(u_end[i]))
            push!(keys, make_var(c.de.equations[j].lhs, τ))
        elseif l isa Number
            push!(u0, l)
            push!(keys, make_var(c.de.equations[j].lhs, τ))
        else
            check = false
            for i=1:length(lhs0)
                l_ = substitute(l, Dict(lhs0[i] => u_end[i]))
                check = !isequal(l_, l)
                check && (push!(u0, l_); push!(keys, make_var(c.de.equations[j].lhs, τ)); break)
            end
            check || error("Could not find initial value for $l !")
        end
    end
    return keys .=> u0
end

"""
    correlation_p0(c::CorrelationFunction, u_end, ps=Pair[])

Find all occurring steady-state values and add them to a list of parameters to
pass this to the `ODEProblem`.

See also: [`CorrelationFunction`](@ref) [`correlation_u0`](@ref)
"""
function correlation_p0(c::CorrelationFunction, u_end, ps=Pair{Any,Any}[])
    if c.steady_state
        steady_vals = c.de0.states
        steady_params = map(_make_parameter, steady_vals)
        ps′ = (ps..., (steady_params .=> u_end)...)
    else
        # Check if <a0> is contained
        avg = average(c.op2_0)
        conj_avg = _conj(avg)
        if avg ∈ Set(c.de0.states)
            p = _make_parameter(avg)
            idx = findfirst(isequal(avg), c.de0.states)
            ps′ = (ps..., p=>u_end[idx])
        elseif conj_avg ∈ Set(c.de0.states)
            p = _make_parameter(conj_avg)
            idx = findfirst(isequal(conj_avg), c.de0.states)
            ps′ = (ps..., p=>u_end[idx])
        else
            ps′ = ps
        end
    end
    return ps′
end

"""
    struct Spectrum

Type representing the spectrum, i.e. the Fourier transform of a
[`CorrelationFunction`](@ref) in steady state.

To actually compute the spectrum at a frequency `ω`, construct the type on top
of a correlation function and call it with `Spectrum(c)(ω,usteady,p0)`.
"""
struct Spectrum
    corr
    Afunc
    bfunc
    cfunc
    A
    b
    c
end

"""
    Spectrum(c::CorrelationFunction, ps=[]; kwargs...)

Create an instance of [`Spectrum`](@ref) corresponding to the Fourier transform
of the [`CorrelationFunction`](@ref) `c`.


Examples
========
```
julia> c = CorrelationFunction(a',a,de;steady_state=true)
⟨a′*a_0⟩

julia> S = Spectrum(c)
ℱ(⟨a′*a_0⟩)(ω)
```
"""
function Spectrum(c::CorrelationFunction, ps=[]; w=Parameter(:ω), kwargs...)
    c.steady_state || error("Cannot use Laplace transform when not in steady state! Use `CorrelationFunction(op1,op2,de0;steady_state=true)` or try computing the Fourier transform of the time evolution of the correlation function directly.")
    de = c.de
    de0 = c.de0
    lhs = getfield.(de.equations, :lhs)
    rhs = getfield.(de.equations, :rhs)
    lhs0 = getfield.(de0.equations, :lhs)
    A,b,c_,Afunc,bfunc,cfunc = _build_spec_func(w, lhs, rhs, c.op2_0, c.op2, lhs0, ps; kwargs...)
    return Spectrum(c, Afunc, bfunc, cfunc, A, b, c_)
end

"""
    (s::Spectrum)(ω::Real,usteady,ps=[];wtol=0)

From an instance of [`Spectrum`](@ref) `s`, actually compute the spectral power
density at the frequency `ω`. Numerically solves the equation `x=inv(A)*b` where
`x` is the vector containing the Fourier transformed correlation function, i.e.
the spectrum is given by `real(x[1])`.
`A` and `b` are a matrix and a vector, respectively, describing the linear system
of equations that needs to be solved to obtain the spectrum.
The tolerance `wtol=0` specifies in which range the frequency should be treated
as zero, i.e. whenever `abs(ω) <= wtol` the term proportional to `1/(im*ω)` is
neglected to avoid divergences.
"""
function (s::Spectrum)(ω::Real,usteady,ps=[];wtol=0)
    A = s.Afunc[1](ω,usteady,ps)
    b = s.bfunc[1](usteady,ps)
    if abs(ω) <= wtol
        b_ = b
    else
        c = s.cfunc[1](ω,usteady,ps)
        b_ = b .+ c
    end
    return 2*real(getindex(A \ b_, 1))
end

"""
    (s::Spectrum)(ω_ls,usteady,ps=[];wtol=0)

From an instance of [`Spectrum`](@ref) `s`, actually compute the spectral power
density at all frequencies in `ω_ls`.
"""
function (s::Spectrum)(ω_ls,usteady,ps=[];wtol=0)
    s_ = Vector{real(eltype(usteady))}(undef, length(ω_ls))
    A = s.Afunc[1](ω_ls[1],usteady,ps)
    b0 = s.bfunc[1](usteady,ps)
    b = copy(b0)
    c = s.cfunc[1](ω_ls[1],usteady,ps)

    if abs(ω_ls[1]) <= wtol
        s_[1] = 2*real(getindex(A \ b, 1))
    else
        s_[1] = 2*real(getindex(A \ (b .+ c), 1))
    end

    Afunc! = (A,ω) -> s.Afunc[2](A,ω,usteady,ps)
    cfunc! = (c,ω) -> s.cfunc[2](c,ω,usteady,ps)
    @inbounds for i=2:length(ω_ls)
        Afunc!(A,ω_ls[i])
        if abs(ω_ls[i]) <= wtol
            s_[i] = 2*real(getindex(A \ b0, 1))
        else
            cfunc!(c,ω_ls[i])
            @. b = b0 + c
            s_[i] = 2*real(getindex(A \ b, 1))
        end
    end

    return s_
end

# Convert to ODESystem
function MTK.ODESystem(c::CorrelationFunction; complete_sys = true, kwargs...)
    τ = MTK.get_iv(c.de)

    ps = []
    for eq∈c.de.equations
        MTK.collect_vars!([],ps,eq.rhs,τ)
    end
    unique!(ps)

    if c.steady_state
        steady_vals = c.de0.states
        steady_hashes = map(hash, steady_vals)
        avg = average(c.op2_0)
        h = hash(avg)
        avg_adj = _adjoint(avg)
        h′ = hash(avg_adj)
        idx = findfirst(isequal(h), steady_hashes)
        if idx === nothing
            idx_ = findfirst(isequal(h′), steady_hashes)
            if idx_ === nothing
                de = c.de
            else
                subs = Dict(average(c.op2) => _adjoint(steady_vals[idx_]))
                de = substitute(c.de, subs)
            end
        else
            subs = Dict(average(c.op2) => steady_vals[idx])
            de = substitute(c.de, subs)
        end
        ps_ = [ps..., steady_vals...]
    else
        avg = average(c.op2_0)
        if avg ∈ Set(c.de0.states)
            avg2 = average(c.op2)
            ps_ = [ps..., avg]
            de = substitute(c.de, Dict(avg2 => avg))
        elseif _conj(avg) ∈ Set(c.de0.states)
            avg2 = average(c.op2)
            ps_ = [ps..., _conj(avg)]
            de = substitute(c.de, Dict(avg2 => avg))
        else
            ps_ = [ps...]
            de = c.de
        end
    end

    ps_avg = filter(x->x isa Average, ps_)
    ps_adj = map(_conj, ps_avg)
    filter!(x->!(x ∈ Set(ps_avg)), ps_adj)
    ps_adj_hash = hash.(ps_adj)

    de_ = deepcopy(de)
    for i=1:length(de.equations)
        lhs = de_.equations[i].lhs
        rhs = substitute_conj(de_.equations[i].rhs, ps_adj, ps_adj_hash)
        de_.equations[i] = Symbolics.Equation(lhs, rhs)
    end

    avg0 = average(c.op2_0)
    if c.steady_state
        steady_params = map(_make_parameter, steady_vals)
        subs_params = Dict(steady_vals .=> steady_params)
        de_ = substitute(de_, subs_params)
    elseif avg0 ∈ Set(c.de0.states)
        avg0_par = _make_parameter(avg0)
        de_ = substitute(de_, Dict(avg0 => avg0_par))
    elseif _conj(avg0) ∈ Set(c.de0.states)
        avg0_par = _make_parameter(_conj(avg0))
        de_ = substitute(de_, Dict(_conj(avg0) => avg0_par))
    end

    eqs = MTK.equations(de_)
    sys = MTK.ODESystem(eqs, τ; kwargs...)
    return complete_sys ? complete(sys) : sys
end

substitute(c::CorrelationFunction, args...; kwargs...) =
    CorrelationFunction(c.op1, c.op2, substitute(c.de0, args...; kwargs...), substitute(c.de, args...; kwargs...))

function _make_parameter(s::Average)
    name = Symbol(string(s))
    sym = Parameter(name)
    return SymbolicUtils.setmetadata(sym, Symbolics.VariableSource, (:_make_parameter, name))
end
_make_parameter(s::SymbolicUtils.Symbolic{<:Parameter}) = s

### Auxiliary functions for CorrelationFunction
function _new_hilbert(h::ProductSpace, aon)
    if length(aon)==1
        return _new_hilbert(h.spaces[aon[1]], 0)
    else
        spaces = [_new_hilbert(h_, 0) for h_ in h.spaces[aon...]]
        return ProductSpace(spaces)
    end
end
_new_hilbert(h::FockSpace, aon) = FockSpace(Symbol(h.name, 0))
_new_hilbert(h::NLevelSpace, aon) = NLevelSpace(Symbol(h.name, 0), h.levels, h.GS)

function _new_operator(op::Destroy, h, aon=op.aon; add_subscript=nothing)
    if isnothing(add_subscript)
        Destroy(h, op.name, aon; op.metadata)
    else
        Destroy(h, Symbol(op.name, :_, add_subscript), aon; op.metadata)
    end
end
function _new_operator(op::Create, h, aon=op.aon; add_subscript=nothing)
    if isnothing(add_subscript)
        Create(h, op.name, aon; op.metadata)
    else
        Create(h, Symbol(op.name, :_, add_subscript), aon; op.metadata)
    end
end
function _new_operator(t::Transition, h, aon=t.aon; add_subscript=nothing)
    if isnothing(add_subscript)
        Transition(h, t.name, t.i, t.j, aon; t.metadata)
    else
        Transition(h, Symbol(t.name, :_, add_subscript), t.i, t.j, aon; t.metadata)
    end
end
_new_operator(x::Number, h, aon=nothing; kwargs...) = x
function _new_operator(t, h, aon=nothing; kwargs...)
    if SymbolicUtils.iscall(t)
        args = []
        if isnothing(aon)
            for arg in SymbolicUtils.arguments(t)
                push!(args, _new_operator(arg, h; kwargs...))
            end
        else
            for arg in SymbolicUtils.arguments(t)
                push!(args, _new_operator(arg,h,aon; kwargs...))
            end
        end
        f = SymbolicUtils.operation(t)
        return f(args...)
    else
        return t
    end
end
function _new_operator(avg::Average, h, aon=nothing; kwargs...)
    op = SymbolicUtils.arguments(avg)[1]
    if isnothing(aon)
        _average(_new_operator(op, h; kwargs...))
    else
        _average(_new_operator(op, h, aon; kwargs...))
    end
end

function _complete_corr!(de,aon0,lhs_new,order,steady_state;
                                mix_choice=maximum,
                                simplify=true,
                                filter_func=nothing,
                                kwargs...)
    vs = de.states
    H = de.hamiltonian
    J = de.jumps
    Jd = de.jumps_dagger
    rates = de.rates

    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)

    vhash_new = map(hash, lhs_new)
    vhash_new′ = map(hash, _adjoint.(lhs_new))
    filter!(!in(vhash_new), vhash_new′)

    function _filter_aon(x) # Filter values that act only on Hilbert space representing system at time t0
        aon = acts_on(x)
        if aon0 in aon
            length(aon)==1 && return false
            return true
        end
        if steady_state # Include terms without t0-dependence only if the system is not in steady state
            h = hash(x)
            return !(h∈vhash_new || h∈vhash_new′)
        else
            return true
        end
    end
    filter!(_filter_aon, missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    filter_reds = de isa ScaledMeanfieldEquations
    filter_reds && filter_redundants!(missed, de.scale_aons, de.names)

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = meanfield(ops_,H,J;
                                Jdagger=Jd,
                                rates=rates,
                                simplify=simplify,
                                order=order,
                                iv=de.iv,
                                kwargs...)

        _append!(de, me)

        vhash_ = hash.(me.states)
        vs′hash_ = hash.(_conj.(me.states))
        append!(vhash, vhash_)
        for i=1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end

        missed = find_missing(me.equations, vhash, vs′hash; get_adjoints=false)
        filter!(_filter_aon, missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
        filter_reds && filter_redundants!(missed, de.scale_aons, de.names)
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        filter!(!filter_func, missed)
        filter_reds && filter_redundants!(missed, de.scale_aons, de.names)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i=1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs)
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end

### Auxiliary functions for Spectrum

function _build_spec_func(ω, lhs, rhs, a1, a0, steady_vals, ps=[])
    s = Dict(a0=>a1)
    ops = [SymbolicUtils.arguments(l)[1] for l in lhs]

    b = [inorder!(average(substitute(op, s))) for op in ops] # Initial values
    c = [SymbolicUtils.simplify(c_ / (1.0im*ω)) for c_ in _find_independent(rhs, a0)]
    aon0 = acts_on(a0)
    @assert length(aon0)==1
    rhs_ = _find_dependent(rhs, aon0[1])
    Ax = [im*ω*lhs[i] - rhs_[i] for i=1:length(lhs)] # Element-wise form of A*x

    # Substitute <a0> by steady-state average <a>
    s_avg = Dict(average(a0) => average(a1))
    Ax = [inorder!(substitute(A, s_avg)) for A∈Ax]
    c = [inorder!(substitute(c_, s_avg)) for c_∈c]

    # Compute symbolic A column-wise by substituting unit vectors into element-wise form of A*x
    A = Matrix{Any}(undef, length(Ax), length(Ax))
    for i=1:length(Ax)
        subs_vals = zeros(length(Ax))
        subs_vals[i] = 1
        subs = Dict(lhs .=> subs_vals)
        A_i = [inorder!(SymbolicUtils.simplify(substitute(Ax[j],subs))) for j=1:length(Ax)]
        A[:,i] = A_i
    end

    # Substitute conjugates
    vs_adj = map(_conj, steady_vals)
    filter!(x->!(x ∈ Set(steady_vals)), vs_adj)
    vs′hash = map(hash, vs_adj)
    A = [substitute_conj(A_,vs_adj,vs′hash) for A_∈A]
    b = [substitute_conj(b_,vs_adj,vs′hash) for b_∈b]
    c = [substitute_conj(c_,vs_adj,vs′hash) for c_∈c]

    # Keep Symbolics.unflatten_long_ops from stepping into symbolic average
    # functions by substituting. This can be removed once averages store
    # operators as metadata
    A1 = [_substitute_vars(A_) for A_∈A]
    b1 = [_substitute_vars(b_) for b_∈b]
    c1 = [_substitute_vars(c_) for c_∈c]
    steady_vals1 = [_substitute_vars(s_) for s_∈steady_vals]

    # Build functions
    Afunc = Symbolics.build_function(A1, ω, steady_vals1, ps; expression=false)
    bfunc = Symbolics.build_function(b1, steady_vals1, ps; expression=false)
    cfunc = Symbolics.build_function(c1, ω, steady_vals1, ps; expression=false)

    return A, b, c, Afunc, bfunc, cfunc
end

function _substitute_vars(t::T) where T <: SymbolicUtils.Symbolic
    if SymbolicUtils.iscall(t)
        f = SymbolicUtils.operation(t)
        if f === sym_average
            sym = Symbol(string(t))
            return SymbolicUtils.setmetadata(SymbolicUtils.Sym{Complex{Real}}(sym),
                Symbolics.VariableSource, (:_substitute_vars, sym))
        else
            args = [_substitute_vars(arg) for arg∈SymbolicUtils.arguments(t)]
            return SymbolicUtils.maketerm(T, f, args, TermInterface.metadata(t))
        end
    else
        return t
    end
end
_substitute_vars(x::Number) = x

_find_independent(rhs::Vector, a0) = [_find_independent(r, a0) for r in rhs]
function _find_independent(r, a0)
    if SymbolicUtils.is_operation(+)(r)
        args_ind = []
        aon0 = acts_on(a0)
        for arg in SymbolicUtils.arguments(r)
            aon = acts_on(arg)
            (aon0 in acts_on(arg) && length(aon)>1) || push!(args_ind, arg)
        end
        isempty(args_ind) && return 0
        return +(args_ind...)
    else
        aon_r = acts_on(r)
        aon0 = acts_on(a0)
        same = (length(aon_r)==length(aon0)) && all(a ∈ aon0 for a∈aon_r) &&
                all(a ∈ aon_r for a ∈ aon0)
        same && return 0
        return r
    end
end

_find_dependent(rhs::Vector, a0) = [_find_dependent(r, a0) for r in rhs]
function _find_dependent(r, a0)
    if SymbolicUtils.is_operation(+)(r)
        args = []
        for arg in SymbolicUtils.arguments(r)
            aon = acts_on(arg)
            (a0 in aon) && length(aon)>1 && push!(args, arg)
        end
        isempty(args) && return 0
        return +(args...)
    else
        aon_r = acts_on(r)
        aon0 = acts_on(a0)
        same = (length(aon_r)==length(aon0)) && all(a ∈ aon0 for a∈aon_r) &&
                all(a ∈ aon_r for a ∈ aon0)
        same || return 0
        return r
    end
end
