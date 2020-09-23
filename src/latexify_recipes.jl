using Latexify
import MacroTools
using LaTeXStrings

const tsym_latex = Ref(:t)
const transition_idx_script = Ref(:^)

"""
    transition_superscript(::Bool)

Specify whether the indices in a [`Transition`](@ref) operator should be
printed as superscript. Default is `true`. If set to `false`, the indices
corresponding to the levels are printed as subscript.
"""
function transition_superscript(x::Bool)
    if x
        transition_idx_script[] = :^
    else
        transition_idx_script[] = :_
    end
    return x
end

function _postwalk_func(x)
    if x==:𝟙
        return "\\mathbb{1}"
    elseif x==:im
        return :i
    elseif MacroTools.@capture(x, dagger(arg_))
        s = "$(arg)^\\dagger"
        return s
    elseif MacroTools.@capture(x, Transition(arg_,i_,j_))
        s = "$(arg)$(transition_idx_script[]){$(i)$(j)}"
        return s
    elseif MacroTools.@capture(x, Sum(arg_, idx_))
        s_arg = latexify(arg; cdot=false)[2:end-1]
        names = getfield.(idx, :name)
        s_names = "$(names[1])"
        for j=2:length(names)
            s_names *= ", $(names[j])"
        end
        s = "\\sum_{$s_names}$(s_arg)"
        return s
    elseif MacroTools.@capture(x, Indexed(arg_, i_))
        s = "$(arg)_{$i}"
        return s
    elseif MacroTools.@capture(x, IndexedParameter(arg_, i_))
        s = "$(arg)_{$(i...)}"
    else
        return x
    end
end

function _postwalk_average(x)
    if MacroTools.@capture(x, AVERAGE(arg_))
        arg = MacroTools.postwalk(_postwalk_func, arg)
        # TODO: clean up; tricky because of nested string conversion of eg average(dagger(a))
        s = string(arg)
        s = replace(s, "\"" => "")
        s = replace(s, "\\\\" => "\\")
        s = replace(s, "*" => "")
        s = "\\langle " * s * "\\rangle "
        return s
    end
    return x
end

@latexrecipe function f(de::AbstractEquation)
    # Options
    env --> :align
    cdot --> false

    lhs, rhs = _latexify(de.lhs,de.rhs)
    return lhs, rhs
end

function _latexify(lhs_::Vector, rhs_::Vector)

    # Convert eqs to Exprs
    rhs = _to_expression.(rhs_)
    rhs = [MacroTools.postwalk(_postwalk_func, eq) for eq=rhs]
    rhs = [MacroTools.postwalk(_postwalk_average, eq) for eq=rhs]
    lhs = _to_expression.(lhs_)
    lhs = [MacroTools.postwalk(_postwalk_func, eq) for eq=lhs]
    lhs = [MacroTools.postwalk(_postwalk_average, eq) for eq=lhs]

    tsym = tsym_latex[]
    dt = Symbol(:d, tsym)
    ddt = :(d/$(dt))
    lhs = [:($ddt * ($l)) for l=lhs]

    return lhs, rhs
end
_latexify(lhs,rhs) = _latexify([lhs],[rhs])

@latexrecipe function f(op::AbstractOperator)
    # Options
    cdot --> false

    ex = _to_expression(op)
    ex = MacroTools.postwalk(_postwalk_func, ex)
    ex = MacroTools.postwalk(_postwalk_average, ex)
    ex = isa(ex,String) ? LaTeXString(ex) : ex
    return ex
end

@latexrecipe function f(s::SymbolicNumber)
    # Options
    cdot --> false

    ex = _to_expression(s)
    ex = MacroTools.postwalk(_postwalk_func, ex)
    ex = MacroTools.postwalk(_postwalk_average, ex)
    ex = isa(ex,String) ? LaTeXString(ex) : ex
    return ex
end
