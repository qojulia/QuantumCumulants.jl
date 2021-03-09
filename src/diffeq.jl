# Required for calculate_tgrad
MTK.detime_dvs(x::Average) = x

# Needed for build_function -- how hacky is this?
MTK.time_varying_as_func(x::Average, sys) = x

MTK.isparameter(::SymbolicUtils.Sym{<:Parameter}) = true

function MTK.ODESystem(he::HeisenbergEquation; ps=nothing, iv=SymbolicUtils.Sym{Real}(:t), kwargs...)
    he.lhs[1] isa Average || error("Cannot convert operator equations to ODESystem. Use `average` first!")
    missed = ps===nothing ? find_missing(he) : find_missing(he;ps=ps)
    isempty(missed) || error("Cannot convert incomplete system to ODESystem. The following averages are missing: $missed")

    vs_adj = map(_conj, he.lhs)
    filter!(x->!_in(x,he.lhs), vs_adj)
    rhs_ = [substitute_conj(r, vs_adj) for r∈he.rhs]

    D = Symbolics.Differential(iv)

    if ps===nothing
        ps′ = []
        for r∈rhs_
            MTK.collect_vars!([],ps′,r,iv)
        end
        unique!(ps′)
    else
        ps′ = ps
    end

    ps_adj = filter(x->x isa Average, ps′)
    if !isempty(ps_adj)
        ps_adj = map(_conj, ps_adj)
        filter!(x->!_in(x,ps′), ps_adj)
        rhs_ = [substitute_conj(r, ps_adj) for r∈rhs_]
    end

    eqs = [Symbolics.Equation(D(l),r) for (l,r)=zip(he.lhs,rhs_)]

    return MTK.ODESystem(eqs, iv, he.lhs, ps′; kwargs...)
end

function substitute_conj(t,vs_adj)
    if SymbolicUtils.istree(t)
        if t isa Average
            if _in(t, vs_adj)
                t′ = _conj(t)
                return conj(t′)
            else
                return t
            end
        else
            _f = x->substitute_conj(x,vs_adj)
            args = map(_f, SymbolicUtils.arguments(t))
            return SymbolicUtils.similarterm(t, SymbolicUtils.operation(t), args)
        end
    else
        return t
    end
end
