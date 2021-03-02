Base.show(io::IO,h::HilbertSpace) = write(io, "ℋ(", h.name, ")")
function Base.show(io::IO,h::ProductSpace)
    show(io, h.spaces[1])
    for i=2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
end
function Base.show(io::IO,h::ClusterSpace)
    write(io, "$(h.N)x")
    show(io, h.original_space)
end

Base.show(io::IO,x::QSym) = write(io, x.name)
Base.show(io::IO,x::Create) = write(io, string(x.name, "′"))
Base.show(io::IO,x::Transition) = write(io, Symbol(x.name,x.i,x.j))

show_brackets = Ref(true)
function Base.show(io::IO,x::QTerm)
    show_brackets[] && write(io,"(")
    show(io, x.arguments[1])
    for i=2:length(x.arguments)
        show(io, x.f)
        show(io, x.arguments[i])
    end
    show_brackets[] && write(io,")")
end

function SymbolicUtils.show_term(io::IO, t::SymbolicUtils.Term{<:AvgSym})
    write(io, "⟨")
    show_brackets[] = false
    show(io, SymbolicUtils.arguments(t)[1])
    show_brackets[] = true
    write(io, "⟩")
end

function Base.show(io::IO,de::AbstractEquation)
    for i=1:length(de)
        write(io, "∂ₜ(")
        show(io, de.lhs[i])
        write(io, ") = ")
        show(io, de.rhs[i])
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

const T_LATEX = Union{<:QNumber,<:AbstractEquation, <:SymbolicUtils.Symbolic{<:CNumber},
        <:CorrelationFunction,<:Spectrum}
Base.show(io::IO, ::MIME"text/latex", x::T_LATEX) = write(io, latexify(x))
