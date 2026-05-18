@latexrecipe function f(eqs::MeanFieldEquations)
    env --> :align
    cdot --> false
    starred --> false
    return eqs.equations
end

@latexrecipe function f(eqs::NoiseMeanFieldEquations)
    env --> :align
    cdot --> false
    starred --> false
    out = Symbolics.Equation[]
    for (de, ne) in zip(eqs.equations, eqs.noise_equations)
        push!(out, de)
        push!(out, ne)
    end
    return out
end
