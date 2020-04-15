global PPRINT = true
set_pprint(val::Bool) = (global PPRINT=val)

import Base: show

function show(stream::IO,a::BasicOperator)
    write(stream,string(a.label))
end
function show(stream::IO,a::Create)
    write(stream,String(a.label))
    write(stream,"ᵗ")
end
function show(stream::IO,a::Expression)
    write(stream, "(")
    for i=1:length(a.args)-1
        show(stream,a.args[i])
        show(stream,a.f)
    end
    show(stream,a.args[end])
    write(stream, ")")
end
show(stream::IO,f::typeof(⊗)) = write(stream,"⊗")
show(stream::IO,a::Identity) = write(stream,"𝟙")
show(stream::IO,::Zero) = show(stream,0)
function show(stream::IO,a::Transition)
    write(stream,string(a.label))
    write(stream,string(a.i))
    write(stream,string(a.j))
end
function show(stream::IO,a::DontSimplify)
    show(stream,a.args[1])
end

function show(stream::IO,i::Index)
    write(stream,string(i.label))
end
function show(stream::IO,a::IndexedOperator{T,<:SymbolicIndex}) where T
    show(stream,a.operator)
    write(stream, "[")
    show(stream,a.index)
    write(stream, "]")
end
function show(stream::IO,a::IndexedOperator{T,<:Int}) where T
    show(stream,a.operator)
    write(stream, "[$(string(a.index))]")
end

function show(stream::IO,s::SumType)
    write(stream, "Σ_"*"$(string(s.f.index.label))( ")
    show(stream,s.args[1])
    write(stream, " )")
end

function show(stream::IO,de::DifferentialEquation)
    write(stream, "∂ₜ(")
    show(stream, de.lhs)
    write(stream, ") = ")
    show(stream, de.rhs)
end
function show(stream::IO,de::DifferentialEquationSet)
    for i=1:length(de.lhs)
        show(stream,DifferentialEquation(de.lhs[i],de.rhs[i]))
        write(stream,"\n")
    end
end

function show(stream::IO, m::MIME"text/latex", x::AbstractOperator)
    if PPRINT
        show(stream, m, sympify(x))
    else
        show(stream, x)
    end
end
function show(stream::IO, m::MIME"text/latex", x::AbstractEquation)
    if PPRINT
        show(stream, m, sympify(x))
    else
        show(stream, x)
    end
end
function show(stream::IO, m::MIME"text/latex", x::Vector{<:DifferentialEquation})
    if PPRINT
        show(stream, m, sympify.(x))
    else
        show(stream, x)
    end
end

function show(stream::IO, m::MIME"text/latex", i::SymbolicIndex)
    if PPRINT
        show(stream, m, sympify(i))
    else
        show(stream, i)
    end
end
