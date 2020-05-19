using Latexify
using MacroTools
using LaTeXStrings

tsym_latex = Ref(:t)

function _postwalk_func(x)
    if x==:𝟙
        return "\\mathbb{1}"
    elseif x==:im
        return :i
    elseif MacroTools.@capture(x, dagger(arg_))
        s = "$(arg)^\\dagger"
        return s
    elseif MacroTools.@capture(x, Transition(arg_,i_,j_))
        s = "$(arg)^{$(i)$(j)}"
        return s
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

@latexrecipe function f(de::DifferentialEquationSet)
    # Options
    env --> :align
    cdot --> false
    # fmt --> "%.3f"

    # Convert eqs to Exprs
    rhs = _to_expression.(de.rhs)
    rhs = [MacroTools.postwalk(_postwalk_func, eq) for eq=rhs]
    rhs = [MacroTools.postwalk(_postwalk_average, eq) for eq=rhs]
    lhs = _to_expression.(de.lhs)
    lhs = [MacroTools.postwalk(_postwalk_func, eq) for eq=lhs]
    lhs = [MacroTools.postwalk(_postwalk_average, eq) for eq=lhs]

    tsym = tsym_latex[]
    dt = Symbol(:d, tsym)
    ddt = :(d/$(dt))
    lhs = [:($ddt * ($l)) for l=lhs]

    return lhs, rhs
end

@latexrecipe function f(de::DifferentialEquation)
    # Options
    env --> :align
    cdot --> false
    # fmt --> "%.3f"

    # Convert eqs to Exprs
    rhs = _to_expression(de.rhs)
    rhs = MacroTools.postwalk(_postwalk_func, rhs)
    rhs = MacroTools.postwalk(_postwalk_average, rhs)
    lhs = _to_expression(de.lhs)
    lhs = MacroTools.postwalk(_postwalk_func, lhs)
    rhs = MacroTools.postwalk(_postwalk_average, lhs)

    tsym = tsym_latex[]
    dt = Symbol(:d, tsym)
    ddt = :(d/$(dt))
    lhs = [:($ddt * ($l)) for l=lhs]

    return lhs, rhs
end


@latexrecipe function f(op::AbstractOperator)
    # Options
    cdot --> false

    ex = _to_expression(op)
    ex = MacroTools.postwalk(_postwalk_func, ex)
    return ex
end

# @latexrecipe function f(avg::AbstractAverage)
#     # Options
#     cdot --> false
#
#     ex = _to_expression(avg)
#     ex = MacroTools.postwalk(_postwalk_func, ex)
#     ex = MacroTools.postwalk(_postwalk_average, ex)
#     return ex
# end

_to_expression(x::Number) = x
function _to_expression(x::Complex)
    iszero(x) && return x
    if iszero(real(x))
        return :( $(imag(x))*im )
    elseif iszero(imag(x))
        return real(x)
    else
        return :( $(real(x)) + $(imag(x))*im )
    end
end
_to_expression(op::BasicOperator) = op.name
_to_expression(op::EmbeddedOperator) = _to_expression(op.operator)
_to_expression(op::Create) = :(dagger($(op.name)))
_to_expression(op::Transition) = :(Transition($(op.name),$(op.i),$(op.j)) )
_to_expression(t::OperatorTerm) = :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )

# function _to_expression(avg::Average)
#     ex = _to_expression(avg.operator)
#     # ex = MacroTools.postwalk(_postwalk_func, ex)
#     return :(AVERAGE($ex))
# end
# _to_expression(t::AverageTerm) = :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
