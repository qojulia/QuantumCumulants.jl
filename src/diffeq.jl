# Relevant parts of ODESystem interface
MTK.independent_variable(he::HeisenbergEquation) = he.iv

MTK.states(he::HeisenbergEquation) = getindex.(he.varmap, 2)

function MTK.equations(he::HeisenbergEquation)
    # Get the MTK variables
    varmap = he.varmap
    vs = []
    for l∈he.lhs
        idx = findfirst(x->isequal(x[1],l),varmap)
        push!(vs, varmap[idx][2])
    end

    # Substitute conjugate variables by explicit conj
    vs_adj = map(_conj, he.lhs)
    filter!(x->!_in(x,he.lhs), vs_adj)
    rhs = [substitute_conj(r, vs_adj) for r∈he.rhs]

    # Substitute to MTK variables on rhs
    subs = Dict(varmap)
    rhs = [substitute(r, subs) for r∈rhs]

    # Return equations
    t = MTK.independent_variable(he)
    D = MTK.Differential(t)
    return [Symbolics.Equation(D(vs[i]), rhs[i]) for i=1:length(vs)]
end

# Substitute conjugate variables
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

# Conversion to ODESystem
MTK.isparameter(::SymbolicUtils.Sym{<:Parameter}) = true

function MTK.ODESystem(he::HeisenbergEquation; kwargs...)
    eqs = MTK.equations(he)
    t = MTK.independent_variable(he)

    # if ps===nothing
    #     ps′ = []
    #     for eq∈eqs
    #         MTK.collect_vars!([],ps′,eq.rhs,t)
    #     end
    #     unique!(ps′)
    # else
    #     ps′ = ps
    # end

    # ps_adj = filter(x->x isa Average, ps′)
    # if !isempty(ps_adj)
    #     ps_adj = map(_conj, ps_adj)
    #     filter!(x->!_in(x,ps′), ps_adj)
    #     eqs = [Symbolics.Equation(eq.lhs, substitute_conj(eq.rhs, ps_adj)) for eq∈eqs]
    # end

    return MTK.ODESystem(eqs; kwargs...)
end
