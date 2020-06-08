Base.show(io::IO,h::HilbertSpace) = write(io, "ℋ(", h.name, ")")
function Base.show(io::IO,h::ProductSpace)
    show(io, h.spaces[1])
    for i=2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
end

function Base.show(io::IO,x::BasicOperator)
    write(io, x.name)
    show_index(io, x.index)
end
function Base.show(io::IO,x::Create)
    write(io, string(x.name, "′"))
    show_index(io, x.index)
end
function Base.show(io::IO,x::Transition)
    write(io, Symbol(x.name,x.i,x.j))
    show_index(io, x.index)
end

show_brackets = Ref(true)
function Base.show(io::IO,x::Union{OperatorTerm,NumberTerm})
    show_brackets[] && write(io,"(")
    show(io, x.arguments[1])
    for i=2:length(x.arguments)
        show(io, x.f)
        show(io, x.arguments[i])
    end
    show_brackets[] && write(io,")")
end

function show_index(io::IO, index)
    if index != default_index()
        write(io, "[")
        show(io, index[1])
        for i=2:length(index)
            write(io, ",")
            show(io, index[i])
        end
        write(io, "]")
    end
end

function Base.show(io::IO, x::Parameter)
    write(io, x.name)
    show_index(io, x.index)
end

function Base.show(io::IO,x::Average)
    write(io,"⟨")
    show_brackets[] = false
    show(io, x.operator)
    show_brackets[] = true
    write(io,"⟩")
end

Base.show(io::IO, i::Index) = Base.show(io, i.i)
function Base.show(io::IO, s::IndexSet)
    write(io, s.name, "(")
    show(io, s.lower:s.upper)
    write(io, ")")
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
