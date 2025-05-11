function average end

"""
    AvgSym <: CNumber

Symbolic number representing the average over an operator.
See also: [`average`](@ref)
"""
struct AvgSym <: CNumber end

const Average = SymbolicUtils.BasicSymbolic{<:AvgSym}

const sym_average = begin # Symbolic function for averages
    T = SymbolicUtils.FnType{Tuple{QNumber}, AvgSym}
    SymbolicUtils.Sym{T}(:avg)
end

# Type promotion -- average(::QNumber)::Number
SymbolicUtils.promote_symtype(::typeof(sym_average), ::Type{<:QNumber}) = AvgSym

# needs a specific symtype overload, otherwise we build the wrong expressions with maketerm
SymbolicUtils.symtype(::T) where T <: Average = QuantumCumulants.AvgSym

# Direct construction of average symbolic expression
function _average(operator)
    return SymbolicUtils.Term{AvgSym}(sym_average, [operator])
end
# ensure that BasicSymbolic{<:AvgSym} are only single averages
function *(a::Average,b::Average)
    if isequal(a,b)
        return SymbolicUtils.Mul(CNumber,1,Dict(a=>2))
    end
    return SymbolicUtils.Mul(CNumber,1,Dict(a=>1,b=>1))
end
function +(a::Average,b::Average)
    if isequal(a,b)
        return SymbolicUtils.Add(CNumber,0,Dict(a=>2))
    end
    return SymbolicUtils.Add(CNumber,0,Dict(a=>1,b=>1))
end

function QuantumAlgebra.acts_on(s::SymbolicUtils.Symbolic)
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        if f === sym_average
            return acts_on(SymbolicUtils.arguments(s)[1])
        else
            aon = []
            for arg in SymbolicUtils.arguments(s)
                append!(aon, acts_on(arg))
            end
            unique!(aon)
            sort!(aon)
            return aon
        end
    else
        return Int[]
    end
end

"""
    average(::QNumber)
    average(::QNumber,order)

Compute the average of an operator. If `order` is given, the [`cumulant_expansion`](@ref)
up to that order is computed immediately.
"""
average(op::QSym) = _average(op)
function average(op::QTerm)
    f = SymbolicUtils.operation(op)
    if f===(Base.:+) || f===(Base.:-) # linearity
        args = map(average, SymbolicUtils.arguments(op))
        return f(args...)
    elseif f === (Base.:*)
        # Move constants out of average
        c = op.arg_c
        op_ = QMul(1,op.args_nc)
        return c*_average(op_)
    else
        error("Unknown function $f")
    end
end

average(x::SNuN) = x
average(x,order;kwargs...) = cumulant_expansion(average(x),order;kwargs...)

function undo_average(t)
    if SymbolicUtils.iscall(t)
        f = SymbolicUtils.operation(t)
        if isequal(f,sym_average) # "===" results in false sometimes in Symbolics version > 5
            return SymbolicUtils.arguments(t)[1]
        else
            args = map(undo_average, SymbolicUtils.arguments(t))
            return f(args...)
        end
    else
        return t
    end
end

function undo_average(eq::Symbolics.Equation)
    lhs = undo_average(eq.lhs)
    rhs = undo_average(eq.rhs)
    return Symbolics.Equation(lhs,rhs)
end

# Standard simplify and expand functions
function SymbolicUtils.simplify(x::QNumber;kwargs...)
    avg = average(x)
    avg_ = SymbolicUtils.simplify(avg;kwargs...)
    return undo_average(avg_)
end

function Symbolics.expand(x::QNumber;kwargs...)
    expansion = average(x)
    expansion_ = SymbolicUtils.expand(expansion; kwargs...)
    return undo_average(expansion_)
end
