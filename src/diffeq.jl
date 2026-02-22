const AbstractIndexedMeanfieldEquations =
    Union{IndexedMeanfieldEquations,EvaledMeanfieldEquations,IndexedMeanfieldNoiseEquations}
# TODO?: create abstract type AbstractIndexedMeanfieldEquations <: AbstractMeanfieldEquations
# which is supertype of all the IndexedEquations

# Relevant parts of System interface
MTK.get_iv(me::AbstractMeanfieldEquations) = me.iv
MTK.unknowns(me::AbstractMeanfieldEquations) = me.states

function MTK.equations(
    me::Union{AbstractMeanfieldEquations,AbstractIndexedMeanfieldEquations},
)
    # Get the MTK variables
    varmap = me.varmap
    vs = MTK.unknowns(me)
    vhash = map(hash, vs)

    # Substitute conjugate variables by explicit conj
    vs′ = map(_conj, vs)
    # vs′ = map(_inconj, vs) # I am not sure why this does not work. 
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
    if me isa AbstractIndexedMeanfieldEquations
        rhs = [substitute_conj_ind(eq.rhs, vs′, vs′hash) for eq ∈ me.equations]
    else
        # For non-indexed equations, we use the standard conjugate substitution
        rhs = [substitute_conj(eq.rhs, vs′, vs′hash) for eq ∈ me.equations]
    end

    # Substitute to MTK variables on rhs
    subs = Dict(varmap)
    rhs = [substitute(r, subs) for r ∈ rhs]
    vs_mtk = getindex.(varmap, 2)

    # Return equations
    t = MTK.get_iv(me)
    D = MTK.Differential(t)
    return [Symbolics.Equation(D(vs_mtk[i]), rhs[i]) for i = 1:length(vs)]
end

# Substitute conjugate variables
function substitute_conj(t::T, vs′, vs′hash) where {T}
    if SymbolicUtils.iscall(t)
        if t isa Average
            if hash(t)∈vs′hash
                t′ = _conj(t)
                return conj(t′)
            else
                return t
            end
        else
            _f = x->substitute_conj(x, vs′, vs′hash)
            args = map(_f, SymbolicUtils.arguments(t))
            return TermInterface.maketerm(
                T,
                SymbolicUtils.operation(t),
                args,
                TermInterface.metadata(t),
            )
        end
    else
        return t
    end
end
function MTK.System(
    me::Union{AbstractMeanfieldEquations,AbstractIndexedMeanfieldEquations},
    iv = me.iv,
    vars = map(last, me.varmap),
    pars = nothing;
    complete_sys = true,
    kwargs...,
)
    eqs = MTK.equations(me)
    pars = isnothing(pars) ? extract_parameters(eqs, iv) : pars
    sys = MTK.System(eqs, iv, vars, pars; kwargs...)
    return complete_sys ? complete(sys) : sys
end


# Substitute conjugate variables for indexed equations
function substitute_conj_ind(t::T, vs′, vs′hash) where {T}
    if SymbolicUtils.iscall(t)
        if t isa Average
            if hash(t)∈vs′hash
                t′ = _inconj(t)
                return conj(t′)
            else
                return t
            end
        else
            _f = x->substitute_conj_ind(x, vs′, vs′hash)
            args = map(_f, SymbolicUtils.arguments(t))
            f = SymbolicUtils.operation(t)
            return TermInterface.maketerm(T, f, args, TermInterface.metadata(t))
        end
    else
        return t
    end
end

function extract_parameters(eqs::Vector{Symbolics.Equation}, iv = nothing)
    params = Set()
    for eq in eqs
        for var in Symbolics.get_variables(eq.rhs)
            if !SymbolicUtils.iscall(var) && !isequal(var, iv)
                push!(params, var)
            end
        end
    end
    return collect(params)
end

function extract_parameters(me::AbstractMeanfieldEquations)
    eqs = MTK.equations(me)

    return extract_parameters(eqs, me.iv)
end
