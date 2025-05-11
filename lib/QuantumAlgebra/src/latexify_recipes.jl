const transition_idx_script = Ref(:^)

"""
    transition_superscript(::Bool)

Specify whether the indices in a [`Transition`](@ref) operator should be
printed as superscript. Default is `true`. If set to `false`, the indices
corresponding to the levels are printed as subscript.
"""
function transition_superscript(x::Bool)
    if x
        transition_idx_script[] = :^
    else
        transition_idx_script[] = :_
    end
    return x
end

function _postwalk_func(x)
    if x==:𝟙
        return "\\mathbb{1}"
    elseif x==:im
        return :i
    elseif MacroTools.@capture(x, dagger(arg_))
        s = "$(arg)^\\dagger"
        return s
    elseif MacroTools.@capture(x, Transition(arg_,i_,j_))
        s = "{$(arg)}$(transition_idx_script[]){{$(i)$(j)}}"
        return s
    elseif MacroTools.@capture(x, Pauli(arg_,i_,))
        s = "{$(arg)}$(transition_idx_script[]){{$(i)}}"
        return s
    elseif MacroTools.@capture(x, Spin(arg_,i_,))
        s = "{$(arg)}$(transition_idx_script[]){{$(i)}}"
        return s
    elseif MacroTools.@capture(x,IndexedVariable(name_,ind_))
        s = "{$name}$(:_){$(ind)}"
        return s
    elseif MacroTools.@capture(x,DoubleIndexedVariable(name_, ind1_, ind2_))
        s = "{$name}$(:_){$(ind1),$(ind2)}"
        return s
    elseif MacroTools.@capture(x, IndexedOperator(op_,in_,i_,j_))
        s = "{$(op)}$(:_){$(in)}$(transition_idx_script[]){{$(i)$(j)}}"
        return s
    elseif MacroTools.@capture(x, IndexedDestroy(op_,in_))
        s = "{$(op)}$(:_){$(in)}"
        return s
    elseif MacroTools.@capture(x, NumberedDestroy(op_,in_))
        s = "{$(op)}$(:_){$(in)}"
        return s
    elseif MacroTools.@capture(x, NumberedOperator(op_,num_,i_,j_))
        s = "{$(op)}$(:_){$(num)}$(transition_idx_script[]){{$(i)$(j)}}"
        return s
    elseif MacroTools.@capture(x,SingleSum(term_, sumInd_, range_, NEI_))
        s = NEI != "" ? "\\underset{$(sumInd)≠$(NEI)}{\\overset{$(range)}{\\sum}} $(_postwalk_func(term))" : "\\underset{$(sumInd)}{\\overset{$(range)}{\\sum}} $(_postwalk_func(term))"
        s = replace(s,"\\\\" => "\\")
        s = replace(s, "\"" => "")
        s = replace(s, "*" => "")
        return s
    elseif MacroTools.@capture(x,DoubleSum(term_, sumInd_, range_, NEI_))
        s = NEI != "" ? "\\underset{$(sumInd)≠$(NEI)}{\\overset{$(range)}{\\sum}} $(_postwalk_func(term))" : "\\underset{$(sumInd)}{\\overset{$(range)}{\\sum}} $(_postwalk_func(term))"
        s = replace(s,"\\\\" => "\\")
        s = replace(s, "\"" => "")
        s = replace(s, "*" => "")
        return s
    elseif MacroTools.@capture(x,IndexedAverageSum(term_, sumInd_, range_, NEI_))
        term = MacroTools.postwalk(_postwalk_average,term)
        s = NEI != "" ? "\\underset{$(sumInd)≠$(NEI)}{\\overset{$(range)}{\\sum}} $(term)" : "\\underset{$(sumInd)}{\\overset{$(range)}{\\sum}} $(term)"
        s = replace(s,"\\\\" => "\\")
        s = replace(s, "\"" => "")
        s = replace(s, "*" => "")
        return s
    elseif MacroTools.@capture(x,IndexedAverageDoubleSum(term_, sumInd_, range_, NEI_))
        term = MacroTools.postwalk(_postwalk_average,term)
        s = NEI != "" ? "\\underset{$(sumInd)≠$(NEI)}{\\overset{$(range)}{\\sum}} $(term)" : "\\underset{$(sumInd)}{\\overset{$(range)}{\\sum}} $(term)"
        s = replace(s,"\\\\" => "\\")
        s = replace(s, "\"" => "")
        s = replace(s, "*" => "")
        return s
    elseif MacroTools.@capture(x,CONJ(arg_))
        arg = MacroTools.postwalk(_postwalk_average,arg)
        s = "{$(arg){^{*}}}"
        return s
    else
        return x
    end
end

function _postwalk_average(x)
    if MacroTools.@capture(x, AVG(arg_))
        arg = MacroTools.postwalk(_postwalk_func, arg)
        # TODO: clean up; tricky because of nested string conversion of eg average(dagger(a))
        s = string(arg)
        s = replace(s, "\"" => "")
        s = replace(s, "\\\\" => "\\")
        s = replace(s, "*" => "")
        s = "\\langle " * s * "\\rangle "
        return s
    end
    return x
end


@latexrecipe function f(op::QNumber)
    # Options
    cdot --> false

    ex = _to_expression(op)
    ex = MacroTools.postwalk(_postwalk_func, ex)
    ex = MacroTools.postwalk(_postwalk_average, ex)
    ex = isa(ex,String) ? latexstring(ex) : ex
    return ex
end

@latexrecipe function f(s::SymbolicUtils.Symbolic{<:CNumber})
    # Options
    cdot --> false

    ex = _to_expression(s)
    ex = MacroTools.postwalk(_postwalk_func, ex)
    ex = MacroTools.postwalk(_postwalk_average, ex)
    ex = isa(ex,String) ? latexstring(ex) : ex
    return ex
end

_to_expression(x::Number) = x
function _to_expression(x::Complex) # For brackets when using latexify
    iszero(x) && return x
    if iszero(real(x))
        return :( $(imag(x))*im )
    elseif iszero(imag(x))
        return real(x)
    else
        return :( $(real(x)) + $(imag(x))*im )
    end
end
function _to_expression(x::Complex{Symbolics.Num}) # forward complex Nums to Symbolics recipes
    iszero(x) && return x
    if iszero(real(x))
        :( $(Symbolics.recipe(imag(x))) * $im )
    elseif iszero(imag(x))
        return :( $(Symbolics.recipe(real(x))) )
    else
        return :( $(Symbolics.recipe(real(x))) + $(Symbolics.recipe(imag(x))) * $im )
    end
end
_to_expression(op::QSym) = op.name
_to_expression(op::Create) = :(dagger($(op.name)))
_to_expression(op::Transition) = :(Transition($(op.name),$(op.i),$(op.j)) )
xyz_sym=[:x,:y,:z]
_to_expression(op::Pauli) = :(Pauli($(op.name),$(xyz_sym[op.axis])))
_to_expression(op::Spin) = :(Spin($(op.name),$(xyz_sym[op.axis])))
function _to_expression(t::QMul)
    args = if SymbolicUtils._isone(t.arg_c)
        t.args_nc
    else
        SymbolicUtils.arguments(t)
    end
    return :( *($(_to_expression.(args)...)) )
end
_to_expression(t::QAdd) = :( +($(_to_expression.(t.arguments)...)) )

_to_expression(p::Parameter) = p.name
function _to_expression(s::SymbolicUtils.Symbolic)
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        fsym = if isequal(f,sym_average) # "===" results in false for symbolics version 5
            :AVG
        elseif f === conj
            :CONJ
        else
            Symbol(f)
        end
        args = map(_to_expression, SymbolicUtils.arguments(s))
        return :( $(fsym)($(args...)) )
    else
        return nameof(s)
    end
end

_to_expression(p::RealParameter) = p.name
