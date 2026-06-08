# Plain-text and LaTeX display for the equation containers, the correlation
# function, and the spectrum. Operator and average rendering is provided by
# SecondQuantizedAlgebra (its `@latexrecipe`s plus the `AvgFunc` hook into
# Symbolics' latexify), so these recipes only assemble the surrounding system.

function Base.show(io::IO, de::AbstractMeanFieldEquations)
    for eq in de.equations
        write(io, "∂ₜ(")
        show(io, eq.lhs)
        write(io, ") = ")
        show(io, eq.rhs)
        write(io, "\n")
    end
    return
end

Base.show(io::IO, c::CorrelationFunction) = show(io, average(c.op1 * c.op2))

function Base.show(io::IO, s::Spectrum)
    write(io, "ℱ(")
    show(io, s.c)
    return write(io, ")(ω)")
end

@latexrecipe function f(de::AbstractMeanFieldEquations)
    env --> :align
    cdot --> false
    D = Symbolics.Differential(de.iv)
    lhs = [D(eq.lhs) for eq in de.equations]
    rhs = [eq.rhs for eq in de.equations]
    return lhs, rhs
end

@latexrecipe function f(de::NoiseMeanFieldEquations)
    env --> :align
    cdot --> false
    D = Symbolics.Differential(de.iv)
    dW, dt = Symbolics.@variables dW dt
    lhs = [D(eq.lhs) for eq in de.equations]
    rhs = [eq.rhs for eq in de.equations] .+ (dW / dt) .* [neq.rhs for neq in de.noise_equations]
    return lhs, rhs
end

@latexrecipe function f(c::CorrelationFunction)
    return average(c.op1 * c.op2)
end

@latexrecipe function f(s::Spectrum)
    inner = latexify(average(s.c.op1 * s.c.op2); env = :raw)
    return latexstring("\\mathcal{F}(", inner, ")(\\omega)")
end

const _LATEX_TYPES = Union{AbstractMeanFieldEquations, CorrelationFunction, Spectrum}
Base.show(io::IO, ::MIME"text/latex", x::_LATEX_TYPES) = write(io, latexify(x))
