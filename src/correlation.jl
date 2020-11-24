struct CorrelationFunction
    op1
    op2
    op2_0
    de0
    de
end

function CorrelationFunction(op1,op2,de0::DifferentialEquation; steady_state=false, add_subscript=0, mix_choice=maximum)
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


    order = maximum(get_order(l) for l in de0.lhs)
    @assert order > 1
    op_ = op1_*op2_
    @assert get_order(op_) <= order

    he = heisenberg(op_,H,J;rates=de0.rates)
    de_ = average(he, order)
    de = _complete_corr(de_, acts_on(op2_), order, steady_state; mix_choice=mix_choice)

    de0_ = DifferentialEquation([_new_operator(l, h) for l in de0.lhs], [_new_operator(r, h) for r in de0.rhs], H, J, de0.rates)
    return CorrelationFunction(op1_, op2_, op2_0, de0_, de)
end

function build_ode(c::CorrelationFunction, ps=[], args...; kwargs...)
    ps_ = (ps..., average(c.op2_0))
    return build_ode(c.de, ps_, args...; kwargs...)
end
generate_ode(c::CorrelationFunction, args...; kwargs...) = Meta.eval(build_ode(c, args...; kwargs...))
substitute(c::CorrelationFunction, args...; kwargs...) =
    CorrelationFunction(c.op1, c.op2, substitute(c.de0, args...; kwargs...), substitute(c.de, args...; kwargs...))

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
        Transition(h, Symbol(t.name, add_subscript), t.i, t.j, aon)
    end
end
_new_operator(x::Number, h, aon=nothing; kwargs...) = x
function _new_operator(t::Union{<:NumberTerm,<:OperatorTerm}, h, aon=nothing; kwargs...)
    args = []
    if isnothing(aon)
        for arg in t.arguments
            push!(args, _new_operator(arg, h; kwargs...))
        end
    else
        for arg in t.arguments
            push!(args, _new_operator(arg,h,aon; kwargs...))
        end
    end
    return t.f(args...)
end
function _new_operator(avg::Average, h, aon=nothing; kwargs...)
    if isnothing(aon)
        Average(_new_operator(avg.operator, h; kwargs...))
    else
        Average(_new_operator(avg.operator, h, aon; kwargs...))
    end
end

function _complete_corr(de,aon0,order,steady_state; mix_choice=maximum)
    lhs = de.lhs
    rhs = de.rhs

    H = de.hamiltonian
    J = de.jumps
    rates = de.rates

    order_lhs = maximum(get_order.(lhs))
    order_rhs = maximum(get_order.(rhs))
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    vs_ = copy(lhs)
    rhs_ = [cumulant_expansion(r, order_) for r in rhs]
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)

    function _filter_aon(x) # Filter values that act only on Hilbert space representing system at time t0
        aon = acts_on(x)
        if aon0 in aon
            length(aon)==1 && return false
            return true
        end
        return !steady_state # Include terms without t0-dependence only if the system is not in steady state
    end
    filter!(_filter_aon, missed)

    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = isempty(J) ? heisenberg(ops,H) : heisenberg(ops,H,J;rates=rates)
        he_avg = average(he,order_;mix_choice=mix_choice)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(x->isa(x,Average),missed)
        filter!(_filter_aon, missed)
    end
    return DifferentialEquation(vs_, rhs_, H, J, rates)
end

struct Spectrum
    corr
    Afunc
    bfunc
end

function Spectrum(c::CorrelationFunction, ps=[]; kwargs...)
    de = c.de
    de0 = c.de0
    Ameta, bmeta = _build_corr_func(de.lhs, de.rhs, c.op2_0, c.op2, de0.lhs, ps; kwargs...)
    Afunc = Meta.eval(Ameta)
    bfunc = Meta.eval(bmeta)
    return Spectrum(c, Afunc, bfunc)
end

function (s::Spectrum)(ω::Real,usteady,ps=[])
    A = s.Afunc(ω,usteady,ps)
    b = s.bfunc(usteady,ps)
    return real(getindex(inv(A)*b, 1))
end

function (s::Spectrum)(ω_ls,usteady,ps=[])
    _Af = ω -> s.Afunc(ω, usteady, ps)
    b = s.bfunc(usteady,ps)
    s_ = Vector{real(eltype(usteady))}(undef, length(ω_ls))
    for i=1:length(ω_ls)
        A = _Af(ω_ls[i])
        s_[i] = real(inv(A)*b)[1]
    end
    return s_
end

function _build_corr_func(lhs, rhs, a1, a0, steady_vals, ps=[]; psym=:p, wsym=:ω, usteady=:usteady)
    s = Dict(a0=>a1)
    ops = getfield.(lhs, :operator)
    lhs_ = [average(substitute(op, s)) for op in ops]

    ω = Parameter{Number}(wsym) # Laplace transform argument i*ω
    b = _find_independent(rhs, a0) # Constant terms; i.e. b in A*x = b
    b = [simplify_constants(b[i] - lhs_[i]) for i=1:length(lhs_)] # Subtract steady-state values
    rhs_ = _find_dependent(rhs, a0)
    Ax = [simplify_constants(- im*ω*lhs[i] + rhs_[i]) for i=1:length(lhs)] # Element-wise form of A*x

    vs = _to_expression.(lhs)
    Ax_ = _to_expression.(Ax)
    b_ = _to_expression.(b)

    # Replace Laplace transform variables
    Ax_ = [MacroTools.postwalk(x -> x in vs ? :( y[$(findfirst(isequal(x), vs))] ) : x, A) for A in Ax_]

    # Replace steady-state values
    if !isempty(steady_vals)
        ss_ = _to_expression.(steady_vals)
        ss_adj = _to_expression.(adjoint.(steady_vals))
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

        # Replace <a0> by <a> steady state value
        a0_ex = _to_expression(average(a0))
        avg1 = average(a1)
        a1_ex = _to_expression(avg1)
        a1_ex_adj = _to_expression(avg1')
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
    end


    # Replace parameters
    if !isempty(ps)
        ps_ = _to_expression.(ps)
        psyms = [:($psym[$i]) for i=1:length(ps)]
        Ax_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, A) for A in Ax_]
        b_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, b1) for b1 in b_]
    end

    line_eqs = [Expr(:(=), :(A[$i,i]), Ax_[i]) for i=1:length(Ax_)]
    ex = build_expr(:block, line_eqs)
    N = length(line_eqs)
    fargs = :($wsym,$usteady,$psym)
    f_ex = :(
        ($fargs) ->
        begin
            T = complex(eltype($usteady))
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

    line_eqs = [Expr(:(=), :(x[$i]), b_[i]) for i=1:length(b_)]
    ex = build_expr(:block, line_eqs)
    N = length(b_)
    fb = :(
        ($usteady,$psym) ->
        begin
            x = zeros(complex(eltype($usteady)), $N)
            begin
                $ex
            end
            return x
        end
    )

    return f_ex, fb
end

_find_independent(rhs::Vector, a0) = [_find_independent(r, a0) for r in rhs]
function _find_independent(r::NumberTerm, a0)
    if r.f === (+)
        args_ind = Number[]
        aon0 = acts_on(a0)
        for arg in r.arguments
            aon0 in acts_on(arg) || push!(args_ind, arg)
        end
        isempty(args_ind) && return 0
        return +(args_ind...)
    else
        return 0
    end
end

_find_dependent(rhs::Vector, a0) = [_find_dependent(r, a0) for r in rhs]
function _find_dependent(r::NumberTerm, a0)
    if r.f === (+)
        args_ind = Number[]
        aon0 = acts_on(a0)
        for arg in r.arguments
            aon0 in acts_on(arg) && push!(args_ind, arg)
        end
        isempty(args_ind) && return 0
        return +(args_ind...)
    else
        return 0
    end
end
