# Relevant parts of ODESystem interface
MTK.get_iv(me::AbstractMeanfieldEquations) = me.iv
MTK.states(me::AbstractMeanfieldEquations) = me.states

function MTK.equations(me::AbstractMeanfieldEquations)
    # Get the MTK variables
    varmap = me.varmap
    vs = MTK.states(me)
    vhash = map(hash, vs)

    # Substitute conjugate variables by explicit conj
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    i = 1
    while i <= length(vs′)
        if vs′hash[i] ∈ vhash
            deleteat!(vs′, i)
            deleteat!(vs′hash, i)
        else
            i += 1
        end
    end
    rhs = [substitute_conj(eq.rhs, vs′, vs′hash) for eq∈me.equations]

    # Substitute to MTK variables on rhs
    subs = Dict(varmap)
    rhs = [substitute(r, subs) for r∈rhs]
    vs_mtk = getindex.(varmap, 2)

    # Return equations
    t = MTK.get_iv(me)
    D = MTK.Differential(t)
    return [Symbolics.Equation(D(vs_mtk[i]), rhs[i]) for i=1:length(vs)]
end

# Substitute conjugate variables
function substitute_conj(t,vs′,vs′hash)
    if SymbolicUtils.istree(t)
        if t isa Average
            if hash(t)∈vs′hash
                t′ = _conj(t)
                return conj(t′)
            else
                return t
            end
        else
            _f = x->substitute_conj(x,vs′,vs′hash)
            args = map(_f, SymbolicUtils.arguments(t))
            return SymbolicUtils.similarterm(t, SymbolicUtils.operation(t), args)
        end
    else
        return t
    end
end

function MTK.ODESystem(me::AbstractMeanfieldEquations, iv=me.iv; kwargs...)
    eqs = MTK.equations(me)
    return MTK.ODESystem(eqs, iv; kwargs...)
end

const AbstractIndexedMeanfieldEquations = Union{IndexedMeanfieldEquations,EvaledMeanfieldEquations}

function MTK.ODESystem(me::AbstractIndexedMeanfieldEquations, iv=me.iv; kwargs...)
    eqs = MTK.equations(me)
    return MTK.ODESystem(eqs, iv; kwargs...)
end

function MTK.equations(me::AbstractIndexedMeanfieldEquations)
    # Get the MTK variables
    varmap = me.varmap
    vs = MTK.states(me)
    vhash = map(hash, vs)

    # Substitute conjugate variables by explicit conj
    vs′ = map(_inconj, vs)
    vs′hash = map(hash, vs′)
    i = 1
    while i <= length(vs′)
        if vs′hash[i] ∈ vhash
            deleteat!(vs′, i)
            deleteat!(vs′hash, i)
        else
            i += 1
        end
    end
    rhs = [substitute_conj_ind(eq.rhs, vs′, vs′hash) for eq∈me.equations]

    # Substitute to MTK variables on rhs
    subs = Dict(varmap)
    rhs = [substitute(r, subs) for r∈rhs]
    vs_mtk = getindex.(varmap, 2)

    # Return equations
    t = MTK.get_iv(me)
    D = MTK.Differential(t)
    return [Symbolics.Equation(D(vs_mtk[i]), rhs[i]) for i=1:length(vs)]
end

# Substitute conjugate variables for indexed equations
function substitute_conj_ind(t,vs′,vs′hash)
    if SymbolicUtils.istree(t)
        if t isa Average
            if hash(t)∈vs′hash
                t′ = _inconj(t)
                return conj(t′)
            else
                return t
            end
        else
            _f = x->substitute_conj_ind(x,vs′,vs′hash)
            args = map(_f, SymbolicUtils.arguments(t))
            return SymbolicUtils.similarterm(t, SymbolicUtils.operation(t), args)
        end
    else
        return t
    end
end
