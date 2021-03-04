"""
    build_function(he::HeisenbergEquation, ps=[]; kwargs...)
"""
function Symbolics.build_function(he::AbstractEquation, ps=[];
                    t=SymbolicUtils.Sym{Real}(:t),
                    iip=true,
                    expression=Val{true},
                    kwargs...)
    rhs, vs = he.rhs, he.lhs
    vs[1] isa Average || error("Cannot build function from q-Number equations!")
    @assert length(rhs) == length(vs)

    vs_adj = map(_conj, vs)
    filter!(x->!_in(x,vs), vs_adj)

    # Check if there are unknown symbols
    missed = find_missing(rhs,vs;vs_adj=vs_adj,ps=ps)
    isempty(missed) || throw_missing_error(missed)

    # Substitute all averages of which the complex conjugate is actually given
    # by an explicit conj
    ps_adj = [_conj(p_) for p_∈ps]
    filter!(x->x isa Average, ps_adj)
    filter!(x->!_in(x,ps), ps_adj)

    rhs_ = [substitute_conj(r,vcat(vs_adj, ps_adj)) for r∈rhs]

    fbuilt = build_function(rhs_, vs, ps, t; expression=expression, kwargs...)
    f = if iip
        fbuilt[2]
    else
        fbuilt[1]
    end
    if expression == Val{true}
        return f
    else
        if iip
            return (du,u,p,t)->f(du,u,p,t)
        else
            return (u,p,t)->f(u,p,t)
        end
    end
end

function throw_missing_error(missed)
    error_msg = "The following parameters or averages are missing: "
    for p1=missed
        error_msg *= "$p1 "
    end
    error(error_msg)
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
