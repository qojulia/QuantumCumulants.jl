import SymPy: sympify

function sympify(a::BasicOperator)
    return SymPy.symbols(string(a.label), commutative=false)
end
function sympify(a::Create)
    return SymPy.symbols(string(a.label), commutative=false)'
end
function sympify(a::Transition)
    return SymPy.symbols(string(a.label)*"^{"*string(a.i)*string(a.j)*"}",commutative=false)
end
sympify(::Identity) = SymPy.Sym(1)#SymPy.symbols("Id",commutative=false)
sympify(::Zero) = SymPy.Sym(0)#SymPy.symbols("Zr",commutative=false)
function sympify(a::NeqIndsProd)
    inds = [a1.index for a1=a.args]
    ind_combs = combinations(inds,2)
    fac = SymPy.Sym(1)
    for ij=ind_combs
        fac *= (1-KroneckerDelta(ij[1],ij[2]))
    end
    return sympify(prod(a.args))*fac
end

sympify(a::Prod) = prod(sympify.(a.args))
function sympify(a::TensorProd)
    args_ = sympify.(a.args)
    return prod(args_)
end
sympify(a::Add) = sum(sympify.(a.args))

function sympify(a::IndexedOperator)
    op_ = sympify(a.operator)
    if isone(op_)
        return 1
    end
    op = SymPy.sympy.IndexedBase(op_)
    i = sympify(a.index)
    return PyCall.py"NCIndexed"(op,i)
end
function sympify(s::SumType)
    arg = sympify(s.args[1])
    i = sympify(s.f.index)
    return SymPy.expand(SymPy.sympy.Sum(arg,(i,i.__pyobject__.lower,i.__pyobject__.upper)))
end

function sympify(de::DifferentialEquation{<:AbstractOperator,<:AbstractOperator};
                tsym=:t)
    d = SymPy.symbols("\\frac{d}{d$(string(tsym))}",commutative=false)
    return SymPy.Eq(d * sympify(de.lhs), sympify(de.rhs))
end
function sympify(de::DifferentialEquation{<:SymPy.Sym,<:SymPy.Sym};
                tsym=:t)
    t = SymPy.symbols(string(tsym))
    lhs_func = SymPy.sympy.Function(string(de.lhs))(t)
    d_lhs = SymPy.sympy.diff(lhs_func, t)
    return SymPy.Eq(d_lhs(t), de.rhs)
end
function sympify(de::DifferentialEquationSet; tsym=:t)
    des = SymPy.Sym[]
    for i=1:length(de.lhs)
        push!(des, sympify(DifferentialEquation(de.lhs[i],de.rhs[i])))
    end
    return des
end

function replace_daggers(s::String)
    if occursin("Dagger",s)
        inds = findfirst("Dagger",s)
        ind0, ind1 = inds[1], inds[end]
        ind2 = findfirst(")",s[ind1:end])[end] + ind1-1
        if ind0>1
            s_ = s[1:ind0-1]
        else
            s_ = ""
        end
        s_ *= s[ind1+1:ind2-1]*"^{\\dagger})"
        if ind2 < length(s)
            s_ *= s[ind2+1:end]
        end
        return replace_daggers(s_)

    elseif occursin("adjoint",s)
        inds = findfirst("adjoint",s)
        ind0, ind1 = inds[1], inds[end]
        ind2 = findfirst(")",s[ind1:end])[end] + ind1-1
        if ind0>1
            s_ = s[1:ind0-1]
        else
            s_ = ""
        end
        s_ *= s[ind1+1:ind2-1]*"^{\\dagger})"
        if ind2 < length(s)
            s_ *= s[ind2+1:end]
        end
        return replace_daggers(s_)
    else
        return replace(s, "*" => "")
    end
end
