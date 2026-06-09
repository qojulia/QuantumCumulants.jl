function Base.show(io::IO, de::AbstractMeanfieldEquations)
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

@latexrecipe function f(de::AbstractMeanfieldEquations)
    env --> :align
    cdot --> false
    D = Symbolics.Differential(de.iv)
    lhs = [D(eq.lhs) for eq in de.equations]
    rhs = [eq.rhs for eq in de.equations]
    return lhs, rhs
end

@latexrecipe function f(de::NoiseMeanfieldEquations)
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

const _LATEX_TYPES = Union{AbstractMeanfieldEquations, CorrelationFunction, Spectrum}
Base.show(io::IO, ::MIME"text/latex", x::_LATEX_TYPES) = write(io, latexify(x))
