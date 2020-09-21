Base.show(io::IO,h::HilbertSpace) = write(io, "ℋ(", h.name, ")")
function Base.show(io::IO,h::ProductSpace)
    show(io, h.spaces[1])
    for i=2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
end

Base.show(io::IO,x::BasicOperator) = write(io, x.name)
Base.show(io::IO,x::Create) = write(io, string(x.name, "′"))
Base.show(io::IO,x::Transition) = write(io, Symbol(x.name,x.i,x.j))

function Base.show(io::IO,x::IndexedDestroy)
    write(io,x.name)
    show(io,x.index)
end
function Base.show(io::IO,x::IndexedCreate)
    write(io,x.name)
    show(io,x.index)
    write(io,"′")
end
Base.show(io::IO,x::IndexedTransition) = write(io, Symbol(x.name,x.i,x.j,x.index.name))
Base.show(io::IO, f::typeof(nip)) = write(io, "⋅")
Base.show(io::IO, x::Index) = write(io, x.name)

function Base.show(io::IO, x::IndexedParameter)
    write(io, x.name)
    write(io, "[")
    show(io, x.index[1])
    for i=2:length(x.index)
        write(io, ", ")
        show(io, x.index[i])
    end
    write(io, "]")
end

show_brackets = Ref(true)
function Base.show(io::IO,x::Union{OperatorTerm,NumberTerm})
    if x.f === Sum
        write(io, "Σ$(x.arguments[2])[")
        show(io, x.arguments[1])
        write(io,"]")
    elseif x.f === (!)
        if x.arguments[1] isa NumberTerm && x.arguments[1].f===(==)
            show(io, x.arguments[1].arguments[1])
            write(io, "≠")
            show(io, x.arguments[1].arguments[2])
        else
            show(io, !)
            write(io, "(")
            for arg in x.arguments
                show(io, arg)
            end
            write(io, ")")
        end
    else
        show_brackets[] && write(io,"(")
        show(io, x.arguments[1])
        for i=2:length(x.arguments)
            show(io, x.f)
            show(io, x.arguments[i])
        end
        show_brackets[] && write(io,")")
    end
end

Base.show(io::IO, x::Parameter) = write(io, x.name)
function Base.show(io::IO,x::Average)
    write(io,"⟨")
    show_brackets[] = false
    show(io, x.operator)
    show_brackets[] = true
    write(io,"⟩")
end

function Base.show(io::IO,de::DifferentialEquation)
    for i=1:length(de)
        write(io, "∂ₜ(")
        show(io, de.lhs[i])
        write(io, ") = ")
        show(io, de.rhs[i])
        write(io, "\n")
    end
end

Base.show(io::IO, ::MIME"text/latex", op::AbstractOperator) = write(io, latexify(op))
Base.show(io::IO, ::MIME"text/latex", de::AbstractEquation) = write(io, latexify(de))
Base.show(io::IO, ::MIME"text/latex", p::SymbolicNumber) = write(io, latexify(p))
