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
function CorrelationFunction(op1,op2,de0::HeisenbergEquation;
                            steady_state=false, add_subscript=0,
                            filter_func=nothing, mix_choice=maximum,
                            iv=SymbolicUtils.Sym{Real}(:τ),
                            simplify=true, kwargs...)
    h1 = hilbert(op1)
    h2 = _new_hilbert(hilbert(op2), acts_on(op2))
    h = h1⊗h2

    H0 = de0.hamiltonian
    J0 = de0.jumps

    op1_ = _new_operator(op1, h)
    op2_ = _new_operator(op2, h, length(h.spaces); add_subscript=add_subscript)
    op2_0 = _new_operator(op2, h)
    H = _new_operator(H0, h)
    J = [_new_operator(j, h) for j in J0]
    lhs_new = [_new_operator(l, h) for l in de0.states]

    order_lhs = maximum(get_order(l) for l in de0.states)
    order_corr = get_order(op1_*op2_)
    order = max(order_lhs, order_corr)
    @assert order > 1
    op_ = op1_*op2_
    @assert get_order(op_) <= order

    de = heisenberg(op_,H,J;rates=de0.rates,iv=iv,expand=true,order=order)
    _complete_corr!(de, length(h.spaces), lhs_new, order, steady_state;
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
        HeisenbergEquation(eqs,eqs_op,lhs_new,ops,H,J,de0.rates,de0.iv,varmap,order)
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
    lhs = [average(substitute(op, subs)) for op in ops]
    u0 = complex(eltype(u_end))[]
    lhs0 = c.de0.states
    τ = MTK.independent_variable(c.de)
    keys = []
    for j=1:length(lhs)
        l=lhs[j]
        l_adj = _adjoint(l)
        if _in(l, lhs0)
            i = findfirst(isequal(l), lhs0)
            push!(u0, u_end[i])
            push!(keys, make_var(c.de.equations[j].lhs, τ))
        elseif _in(l_adj, lhs0)
            i = findfirst(isequal(l_adj), lhs0)
            push!(u0, conj(u_end[i]))
            push!(keys, make_var(c.de.equations[j].lhs, τ))
        else
            check = false
            for i=1:length(lhs0)
                l_ = substitute(l, Dict(lhs0[i] => u_end[i]))
                check = !isequal(l_, l)
                check && (push!(u0, l_); push!(keys, make_var(c.de.equations[i].lhs, τ)); break)
            end
            check || error("Could not find initial value for $l !")
        end
    end
    return keys .=> u0
end

"""
    correlation_p0(c::CorrelationFunction, u_end, ps=[])

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
        if _in(avg, c.de0.states)
            p = _make_parameter(avg)
            idx = findfirst(isequal(avg), c.de0.states)
            ps′ = (ps..., p=>u_end[idx])
        elseif _in(conj_avg, c.de0.states)
            p = _make_parameter(conj_avg)
            idx = findfirst(isequal(conj_avg), c.de0.states)
            ps′ = (ps..., p=>u_end[idx])
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
struct Spectrum{C,FA,FB,FAS,FBS}
    corr::C
    Afunc::FA
    bfunc::FB
    Asym::FAS
    bsym::FBS
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
function Spectrum(c::CorrelationFunction, ps=[]; kwargs...)
    c.steady_state || error("Cannot use Laplace transform when not in steady state! Use `CorrelationFunction(op1,op2,de0;steady_state=true)` or try computing the Fourier transform of the time evolution of the correlation function directly.")
    de = c.de
    de0 = c.de0
    fAsym, fbsym, Ameta, bmeta = _build_spec_func(de.lhs, de.rhs, c.op2_0, c.op2, de0.lhs, ps; kwargs...)
    Afunc = Meta.eval(Ameta)
    bfunc = Meta.eval(bmeta)
    Asymfunc = Meta.eval(fAsym)
    bsymfunc = Meta.eval(fbsym)
    return Spectrum(c, Afunc, bfunc, Asymfunc, bsymfunc)
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
    A = s.Afunc(ω,usteady,ps)
    b = s.bfunc(ω,usteady,ps,wtol)
    return 2*real(getindex(inv(A)*b, 1))
end

"""
    (s::Spectrum)(ω_ls,usteady,ps=[];wtol=0)

From an instance of [`Spectrum`](@ref) `s`, actually compute the spectral power
density at all frequencies in `ω_ls`.
"""
function (s::Spectrum)(ω_ls,usteady,ps=[];wtol=0)
    _Af = ω -> s.Afunc(ω, usteady, ps)
    _bf = ω -> s.bfunc(ω, usteady, ps, wtol)
    s_ = Vector{real(eltype(usteady))}(undef, length(ω_ls))
    for i=1:length(ω_ls)
        A = _Af(ω_ls[i])
        b = _bf(ω_ls[i])
        s_[i] = 2*real(inv(A)*b)[1]
    end
    return s_
end


"""
    (s::Spectrum)(ω::Symbolic,ps=[])
    (s::Spectrum)(ω::Symbolic,steady_vals,ps)

Compute the symbolic system of linear equations that is solved numerically when
computing the spectrum. Returns a symbolic version of the matrix `A` and the
vector `b` describing the linear system of equations that needs to be solved
to obtain the spectrum.
"""
function (s::Spectrum)(ω::SymbolicUtils.Symbolic,ps=[])
    steady_vals = s.corr.de0.lhs
    return s(ω,steady_vals,ps)
end
function (s::Spectrum)(ω::SymbolicUtils.Symbolic,steady_vals,ps)
    A = s.Asym(ω,steady_vals,ps)
    b = s.bsym(ω,steady_vals,ps)
    return qsimplify.(A), qsimplify.(b)
end

# Convert to ODESystem
function MTK.ODESystem(c::CorrelationFunction; kwargs...)
    τ = MTK.independent_variable(c.de)

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
        if _in(avg, c.de0.states) || _in(_conj(avg), c.de0.states)
            ps_ = [ps..., average(c.op2)]
        else
            ps_ = [ps...]
        end
        de = c.de
    end

    ps_avg = filter(x->x isa Average, ps_)
    ps_adj = map(_conj, ps_avg)
    filter!(x->!_in(x,ps_avg), ps_adj)
    ps_adj_hash = hash.(ps_adj)

    de_ = deepcopy(de)
    for i=1:length(de.equations)
        lhs = de_.equations[i].lhs
        rhs = substitute_conj(de_.equations[i].rhs, ps_adj, ps_adj_hash)
        de_.equations[i] = Symbolics.Equation(lhs, rhs)
    end

    if c.steady_state
        steady_params = map(_make_parameter, steady_vals)
        subs_params = Dict(steady_vals .=> steady_params)
        for i=1:length(de.equations)
            de_.equations[i] = substitute(de_.equations[i], subs_params)
        end
    end

    ps_ = map(_make_parameter, ps_)

    eqs = MTK.equations(de_)
    return MTK.ODESystem(eqs, τ; kwargs...)
end

substitute(c::CorrelationFunction, args...; kwargs...) =
    CorrelationFunction(c.op1, c.op2, substitute(c.de0, args...; kwargs...), substitute(c.de, args...; kwargs...))

function _make_parameter(s::Average)
    name = Symbol(string(s))
    return SymbolicUtils.Sym{Parameter}(name)
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
        Destroy(h, op.name, aon)
    else
        Destroy(h, Symbol(op.name, :_, add_subscript), aon)
    end
end
function _new_operator(op::Create, h, aon=op.aon; add_subscript=nothing)
    if isnothing(add_subscript)
        Create(h, op.name, aon)
    else
        Create(h, Symbol(op.name, :_, add_subscript), aon)
    end
end
function _new_operator(t::Transition, h, aon=t.aon; add_subscript=nothing)
    if isnothing(add_subscript)
        Transition(h, t.name, t.i, t.j, aon)
    else
        Transition(h, Symbol(t.name, :_, add_subscript), t.i, t.j, aon)
    end
end
_new_operator(x::Number, h, aon=nothing; kwargs...) = x
function _new_operator(t, h, aon=nothing; kwargs...)
    if SymbolicUtils.istree(t)
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
function _new_operator(avg::SymbolicUtils.Term{<:AvgSym}, h, aon=nothing; kwargs...)
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
    rates = de.rates

    order_lhs = maximum(get_order.(vs))
    order_rhs = 0
    for i=1:length(de.equations)
        k = get_order(de.equations[i].rhs)
        k > order_rhs && (order_rhs = k)
    end
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

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

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        he = heisenberg(ops_,de.hamiltonian,de.jumps;
                                rates=de.rates,
                                simplify=simplify,
                                expand=true,
                                order=order_,
                                iv=de.iv,
                                kwargs...)

        _append!(de, he)

        vhash_ = hash.(he.states)
        vs′hash_ = hash.(_conj.(he.states))
        append!(vhash, vhash_)
        for i=1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end

        missed = find_missing(he.equations, vhash, vs′hash; get_adjoints=false)
        filter!(_filter_aon, missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        filter!(!filter_func, missed)
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

function _build_spec_func(lhs, rhs, a1, a0, steady_vals, ps=[]; psym=:p, wsym=:ω, usteady=:usteady)
    s = Dict(a0=>a1)
    ops = [SymbolicUtils.arguments(l)[1] for l in lhs]

    ω = Parameter(wsym) # Laplace transform argument i*ω
    b = [average(qsimplify(substitute(op, s))) for op in ops] # Initial values
    c = [qsimplify(c_ / (1.0im*ω)) for c_ in _find_independent(rhs, a0)]
    aon0 = acts_on(a0)
    @assert length(aon0)==1
    rhs_ = _find_dependent(rhs, aon0[1])
    Ax = [qsimplify(im*ω*lhs[i] - rhs_[i]) for i=1:length(lhs)] # Element-wise form of A*x

    vs = _to_expression.(lhs)
    Ax_ = _to_expression.(Ax)
    b_ = _to_expression.(b)
    c_ = _to_expression.(c)

    # Replace Laplace transform variables
    Ax_ = [MacroTools.postwalk(x -> x in vs ? :( y[$(findfirst(isequal(x), vs))] ) : x, A) for A in Ax_]

    # Replace steady-state values
    if !isempty(steady_vals)
        ss_ = _to_expression.(steady_vals)
        ss_adj = _to_expression.(_adjoint.(steady_vals))
        ssyms = [:($usteady[$i]) for i=1:length(steady_vals)]
        _pw = function(x)
            if x in ss_
                ssyms[findfirst(isequal(x), ss_)]
            elseif x in ss_adj
                :( conj($(ssyms[findfirst(isequal(x), ss_adj)])) )
            else
                x
            end
        end
        Ax_ = [MacroTools.postwalk(_pw, A) for A in Ax_]
        b_ = [MacroTools.postwalk(_pw, b1) for b1 in b_]
        c_ = [MacroTools.postwalk(_pw, c1) for c1 in c_]

        # Replace <a0> by <a> steady state value
        a0_ex = _to_expression(average(a0))
        avg1 = average(a1)
        a1_ex = _to_expression(avg1)
        a1_ex_adj = _to_expression(_adjoint(avg1))
        _pw2 = function(x)
            if x == a0_ex
                i = findfirst(isequal(a1_ex), ss_)
                if isnothing(i)
                    j = findfirst(isequal(a1_ex_adj), ss_)
                    return :( conj($(ssyms[j])) )
                else
                    return ssyms[i]
                end
            else
                x
            end
        end
        Ax_ = [MacroTools.postwalk(_pw2, A) for A in Ax_]
        b_ = [MacroTools.postwalk(_pw2, b1) for b1 in b_]
        c_ = [MacroTools.postwalk(_pw2, c1) for c1 in c_]
    end


    # Replace parameters
    if !isempty(ps)
        ps_ = _to_expression.(ps)
        psyms = [:($psym[$i]) for i=1:length(ps)]
        Ax_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, A) for A in Ax_]
        b_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, b1) for b1 in b_]
        c_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, c1) for c1 in c_]
    end

    # Obtain A
    line_eqs = [Expr(:(=), :(A[$i,i]), Ax_[i]) for i=1:length(Ax_)]
    ex = Expr(:block, line_eqs...)
    N = length(line_eqs)

    # Function for building numeric A
    fargs = :($wsym,$usteady,$psym)
    fA = :(
        ($fargs) ->
        begin
            T = complex(promote_type(map(eltype, $fargs)...))
            A = Matrix{T}(undef, $N, $N)
            y = zeros(T, $N)
            for i=1:$N
                y[i] = one(T)
                begin
                    $ex
                end
                y[i] = zero(T)
            end
            return A
        end
    )
    # Function for building symbolic A
    fAsym = :(
        ($fargs) ->
        begin
            A = Matrix{Any}(undef, $N, $N)
            y = zeros(Number, $N)
            for i=1:$N
                y[i] = 1
                begin
                    $ex
                end
                y[i] = 0
            end
            return A
        end
    )

    # Obtain b
    line_eqs = [Expr(:(=), :(x[$i]), b_[i]) for i=1:length(b_)]
    ex0 = Expr(:block, line_eqs...)
    eqs_nz = [Expr(:(=), :(x[$i]), :($(b_[i]) + $(c_[i]))) for i=1:length(c_)]
    ex_nz = Expr(:block, eqs_nz...)
    N = length(b_)
    # Function for numeric b
    fb = :(
        ($wsym,$usteady,$psym,wtol=0) ->
        begin
            T = complex(promote_type(eltype($usteady), eltype($psym)))
            x = zeros(T, $N)
            if abs($wsym)<=wtol
                $ex0
            else
                $ex_nz
            end
            return x
        end
    )
    # Function for symbolic b
    fbsym = :(
        ($fargs) ->
        begin
            x = zeros(Number, $N)
            begin
                $ex_nz
            end
            return x
        end
    )

    return fAsym, fbsym, fA, fb
end

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
        return 0
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
        return 0
    end
end
