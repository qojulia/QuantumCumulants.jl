using Latexify
import MacroTools
using LaTeXStrings

@latexrecipe function f(de::AbstractMeanfieldEquations)
    # Options
    env --> :align
    cdot --> false

    lhs = getfield.(de.equations, :lhs)
    rhs = getfield.(de.equations, :rhs)
    lhs, rhs = _latexify(lhs, rhs, de.iv)
    return lhs, rhs
end

@latexrecipe function f(de::Union{MeanfieldNoiseEquations,IndexedMeanfieldNoiseEquations})
    # Options
    env --> :align
    cdot --> false

    lhs = getfield.(de.equations, :lhs)
    rhs = getfield.(de.equations, :rhs)
    rhs_noise = getfield.(de.noise_equations, :rhs)
    npr = cnumbers("dW/dt")
    lhs, rhs = QuantumCumulants._latexify(lhs, rhs.+npr.*rhs_noise, de.iv)
    return lhs, rhs
end

function _latexify(lhs_, rhs_, t)

    # Convert eqs to Exprs
    rhs = _to_expression.(rhs_)
    rhs = [MacroTools.postwalk(SQA._postwalk_func, eq) for eq=rhs]
    rhs = [MacroTools.postwalk(SQA._postwalk_average, eq) for eq=rhs]
    lhs = _to_expression.(lhs_)
    lhs = [MacroTools.postwalk(SQA._postwalk_func, eq) for eq=lhs]
    lhs = [MacroTools.postwalk(SQA._postwalk_average, eq) for eq=lhs]

    tsym = Symbol(t)
    dt = Symbol(:d, tsym)
    ddt = :(d/$(dt))
    lhs = [:($ddt * ($l)) for l=lhs]

    return lhs, rhs
end
_latexify(lhs,rhs) = _latexify([lhs],[rhs])

@latexrecipe function f(c::CorrelationFunction)
    avg = average(c.op1*c.op2)
    return latexify(avg)
end

@latexrecipe function f(S::Spectrum)
    l = latexstring(L"\mathcal{F}(", latexify(S.corr), L")(\omega)")
    return l
end

function SQA._to_expression(s::SymbolicUtils.Symbolic)
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        fsym = if isequal(f, sym_average) # "===" results in false for symbolics version 5
            :AVG
        elseif f === conj
            :CONJ
        else
            Symbol(f)
        end
        args = map(_to_expression, SymbolicUtils.arguments(s))
        return :( $(fsym)($(args...)) )
    else
        return nameof(s)
    end
end
