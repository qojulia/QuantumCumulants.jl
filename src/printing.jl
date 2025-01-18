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
Base.show(io::IO,x::Pauli) = write(io, Symbol(x.name,xyz_sym[x.axis]))
Base.show(io::IO,x::Spin) = write(io, Symbol(x.name,xyz_sym[x.axis]))

show_brackets = Ref(true)
function Base.show(io::IO,x::QTerm)
    show_brackets[] && write(io,"(")
    show(io, SymbolicUtils.arguments(x)[1])
    f = SymbolicUtils.operation(x)
    for i=2:length(SymbolicUtils.arguments(x))
        show(io, f)
        show(io, SymbolicUtils.arguments(x)[i])
    end
    show_brackets[] && write(io,")")
end

function Base.show(io::IO,x::QMul)
    if !SymbolicUtils._isone(x.arg_c)
        show(io, x.arg_c)
        show(io, *)
    end
    show_brackets[] && write(io, "(")
    show(io, x.args_nc[1])
    for i=2:length(x.args_nc)
        show(io, *)
        show(io, x.args_nc[i])
    end
    show_brackets[] && write(io, ")")
end

function SymbolicUtils.show_term(io::IO, t::Average)
    write(io, "⟨")
    show_brackets[] = false
    show(io, SymbolicUtils.arguments(t)[1])
    show_brackets[] = true
    write(io, "⟩")
end

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

const T_LATEX = Union{<:QNumber,<:AbstractMeanfieldEquations, <:SymbolicUtils.Symbolic{<:CNumber},
        <:CorrelationFunction,<:Spectrum}
Base.show(io::IO, ::MIME"text/latex", x::T_LATEX) = write(io, latexify(x))


function Base.show(io::IO, s::Sum)
    write(io, "Σ(")
    show(io, s.term)
    write(io, ", ")
    show(io, s.index)
    write(io, ")")
end

function SymbolicUtils.show_term(io::IO, s::SymbolicUtils.Symbolic{<:CSumSym})
    write(io, "Σ(")
    term, index = SymbolicUtils.arguments(s)
    show(io, term)
    write(io, ", ")
    show(io, index)
    write(io, ")")
end

function SymbolicUtils.show_term(io::IO, p::SymbolicUtils.Symbolic{<:IndexedParameterSym})
    args = SymbolicUtils.arguments(p)
    show(io, args[1])
    write(io, "[")
    for i=2:length(args) - 1
        show(io, args[i])
        write(io, ", ")
    end
    if length(args) > 2
        show(io, args[end])
    end
    write(io, "]")
end


Base.show(io::IO, index::Index) = write(io, index.name)

function Base.show(io::IO,op::IndexedOperator)
    op_ = op.op
    if typeof(op_) <:Transition
        write(io,Symbol(op_.name,op_.i,op_.j,op.ind.name))
    elseif op_ isa Destroy
        write(io,Symbol(op_.name,op.ind.name))
    elseif op_ isa Create
        write(io,Symbol(op_.name,op.ind.name,"'"))
    else
        write(io,op_.name)
    end
end
