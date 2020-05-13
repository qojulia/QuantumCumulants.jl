using Latexify
using MacroTools

tsym_latex = Ref(:t)
simplified_tensor_latex = Ref(true)
@latexrecipe function f(de::DifferentialEquationSet)
    # Options
    env --> :align
    cdot --> false
    # fmt --> "%.3f"

    # Convert eqs to Exprs
    rhs = _to_expression.(de.rhs)
    rhs = [MacroTools.postwalk(x -> MacroTools.@capture(x, dagger(arg_)) ? "$(arg)^\\dagger" : x, eq) for eq=rhs]
    rhs = [MacroTools.postwalk(x -> x==:𝟙 ? "\\mathbb{1}" : x, eq) for eq=rhs]

    lhs = _to_expression.(de.lhs)
    lhs = [MacroTools.postwalk(x -> MacroTools.@capture(x, dagger(arg_)) ? "$(arg)^\\dagger" : x, eq) for eq=lhs]
    lhs = [MacroTools.postwalk(x -> x==:𝟙 ? "\\mathbb{1}" : x, eq) for eq=lhs]
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
    rhs = MacroTools.postwalk(x -> MacroTools.@capture(x, dagger(arg_)) ? "$(arg)^\\dagger" : x, rhs)
    rhs = MacroTools.postwalk(x -> x==:𝟙 ? "\\mathbb{1}" : x, rhs)

    lhs = _to_expression(de.lhs)
    lhs = MacroTools.postwalk(x -> MacroTools.@capture(x, dagger(arg_)) ? "$(arg)^\\dagger" : x, lhs)
    lhs = MacroTools.postwalk(x -> x==:𝟙 ? "\\mathbb{1}" : x, lhs)
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
    ex = MacroTools.postwalk(x -> MacroTools.@capture(x, dagger(arg_)) ? "$(arg)^\\dagger" : x, ex)
    ex = MacroTools.postwalk(x -> x==:𝟙 ? "\\mathbb{1}" : x, ex)
    return ex
end

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
_to_expression(op::Create) = :(dagger($(op.name)))
_to_expression(t::OperatorTerm) = :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
function _to_expression(t::OperatorTerm{<:typeof(⊗)})
    if simplified_tensor_latex[]
        # Replace ⊗ by * and simplify to get rid of identities
        t_ = simplify_operators(OperatorTerm(*, t.arguments))
        return _to_expression(t_)
    else
        return :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
    end
end
