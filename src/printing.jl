using .SecondQuantizedAlgebra: show_brackets

function Base.show(io::IO,de::AbstractMeanfieldEquations)
    for i=1:length(de.equations)
        write(io, "∂ₜ(")
        show(io, de.equations[i].lhs)
        write(io, ") = ")
        show(io, de.equations[i].rhs)
        write(io, "\n")
    end
end

function Base.show(io::IO, c::CorrelationFunction)
    show(io, average(c.op1*c.op2))
end
function Base.show(io::IO, S::Spectrum)
    write(io, "ℱ(")
    show(io, S.corr)
    write(io, ")(ω)")
end

Base.show(io::IO, c::CallableTransition) = write(io, c.name)

const T_LATEX = Union{<:AbstractMeanfieldEquations,<:CorrelationFunction,<:Spectrum}
Base.show(io::IO, ::MIME"text/latex", x::T_LATEX) = write(io, latexify(x))
