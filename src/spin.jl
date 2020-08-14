"""
    SpinSpace <: HilbertSpace

[`HilbertSpace`](@ref) defining a space for spin operators.
"""
struct SpinSpace{S,N} <: HilbertSpace
    name::S
    spin::N
end
Base.:(==)(h1::T,h2::T) where T<:SpinSpace = (h1.name==h2.name && h1.spin==h2.spin)
Base.hash(h::SpinSpace, i::UInt) = hash(h.name, hash(h.spin, i))
isspinN(N,h::SpinSpace,args...) = (N==h.spin)
isspinN(N,h::ProductSpace,aon) = isspinN(N,h.spaces[aon])
isspinN(N,x) = false

abstract type SpinOperator <: BasicOperator end
abstract type SymmetricSpinOperator <: BasicOperator end

isspinop(::SpinOperator) = true
is_symspin(::SymmetricSpinOperator) = true
isanyspin(::Union{SpinOperator,SymmetricSpinOperator}) = true

isspinop(::SymbolicUtils.Sym{<:SpinOperator}) = true
is_symspin(::SymbolicUtils.Sym{<:SymmetricSpinOperator}) = true
isanyspin(::SymbolicUtils.Sym{<:Union{SpinOperator,SymmetricSpinOperator}}) = true


"""
    SigmaX <: SpinOperator

Operator on a [`SpinSpace`](@ref) representing the spin x component.
"""
struct SigmaX{H<:HilbertSpace,S,A} <: SpinOperator
    hilbert::H
    name::S
    aon::A
    function SigmaX{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SigmaX}(gensym(:SigmaX))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issigmax(a::SymbolicUtils.Sym{T}) where {T<:SigmaX} = true

"""
    SigmaY <: SpinOperator

Operator on a [`SpinSpace`](@ref) representing the spin y component.
"""
struct SigmaY{H<:HilbertSpace,S,A} <: SpinOperator
    hilbert::H
    name::S
    aon::A
    function SigmaY{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SigmaY}(gensym(:SigmaY))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issigmay(a::SymbolicUtils.Sym{T}) where {T<:SigmaY} = true

"""
    SigmaZ <: SpinOperator

Operator on a [`SpinSpace`](@ref) representing the spin y component.
"""
struct SigmaZ{H<:HilbertSpace,S,A} <: SpinOperator
    hilbert::H
    name::S
    aon::A
    function SigmaZ{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SigmaZ}(gensym(:SigmaZ))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issigmaz(a::SymbolicUtils.Sym{T}) where {T<:SigmaZ} = true

"""
    SymmetrizedSigmaX <: SymmetricSpinOperator

Operator on a [`SpinSpace`](@ref) representing the spin x component.
"""
struct SymmetrizedSigmaX{H<:HilbertSpace,S,A} <: SymmetricSpinOperator
    hilbert::H
    name::S
    aon::A
    function SymmetrizedSigmaX{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SymmetrizedSigmaX}(gensym(:SymmetrizedSigmaX))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
is_symmetrized_sigmax(a::SymbolicUtils.Sym{T}) where {T<:SymmetrizedSigmaX} = true

"""
    SymmetrizedSigmaY <: SymmetricSpinOperator

Operator on a [`SpinSpace`](@ref) representing the spin y component.
"""
struct SymmetrizedSigmaY{H<:HilbertSpace,S,A} <: SymmetricSpinOperator
    hilbert::H
    name::S
    aon::A
    function SymmetrizedSigmaY{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SymmetrizedSigmaY}(gensym(:SymmetrizedSigmaY))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
is_symmetrized_sigmay(a::SymbolicUtils.Sym{T}) where {T<:SymmetrizedSigmaY} = true

"""
    SymmetrizedSigmaZ <: SymmetricSpinOperator

Operator on a [`SpinSpace`](@ref) representing the spin y component.
"""
struct SymmetrizedSigmaZ{H<:HilbertSpace,S,A} <: SymmetricSpinOperator
    hilbert::H
    name::S
    aon::A
    function SymmetrizedSigmaZ{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SymmetrizedSigmaZ}(gensym(:SymmetrizedSigmaZ))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
is_symmetrized_sigmaz(a::SymbolicUtils.Sym{T}) where {T<:SymmetrizedSigmaZ} = true

symmetrize(x) = x
unsymmetrize(x) = x
for f in [:SigmaX,:SigmaY,:SigmaZ]
    fsym = Symbol(:Symmetrized,f)
    @eval symmetrize(op::($f)) = $(fsym)(op.hilbert,op.name,op.aon)
    @eval unsymmetrize(op::($fsym)) = $(f)(op.hilbert,op.name,op.aon)
end
for f in [:symmetrize,:unsymmetrize]
    @eval $(f)(s::SymbolicUtils.Symbolic) = _to_symbolic(($f)(_to_qumulants(s)))
    @eval function $(f)(x::OperatorTerm)
        args = [unsymmetrize(arg) for arg in x.arguments]
        return x.f(args...)
    end
end

for f in [:SigmaX,:SigmaY,:SigmaZ,:SymmetrizedSigmaX,:SymmetrizedSigmaY,:SymmetrizedSigmaZ]
    @eval $(f)(hilbert::H,name::S,aon::A) where {H,S,A} = $(f){H,S,A}(hilbert,name,aon)
    @eval $(f)(hilbert::SpinSpace,name) = $(f)(hilbert,name,1)
    @eval function $(f)(hilbert::ProductSpace,name)
        i = findall(x->isa(x,SpinSpace),hilbert.spaces)
        if length(i)==1
            return $(f)(hilbert,name,i[1])
        else
            isempty(i) && error("Can only create $($(f)) on SpinSpace! Not included in $(hilbert)")
            length(i)>1 && error("More than one SpinSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!")
        end
    end
    @eval function embed(h::ProductSpace,op::T,aon::Int) where T<:($(f))
        check_hilbert(h.spaces[aon],op.hilbert)
        op_ = $(f)(h,op.name,aon)
        return op_
    end
    @eval function Base.hash(op::T, h::UInt) where T<:($(f))
        hash(op.hilbert, hash(op.name, hash(op.aon, hash($(f), h))))
    end
    @eval Base.adjoint(op::T) where T<:($(f)) = op
    @eval isspinN(N,op::T) where T<:($(f)) = isspinN(N,op.hilbert,op.aon)
end

# Commutation relation in simplification
for (S,Ssym,l) in zip(permutations([:SigmaX,:SigmaY,:SigmaZ]),permutations([:SymmetrizedSigmaX,:SymmetrizedSigmaY,:SymmetrizedSigmaZ]), permutations([1,2,3]))
    ϵ = levicivita(l)
    @eval function commute_spin(x::SymbolicUtils.Sym{S1}, y::SymbolicUtils.Sym{S2}) where {S1<:($(S[1])),S2<:($(S[2]))}
        op_x = _to_qumulants(x)
        z = _to_symbolic($(Ssym[3])(op_x.hilbert,op_x.name,op_x.aon))
        x_ = symmetrize(x)
        y_ = symmetrize(y)
        return x_*y_ + y_*x_ + (im*$(ϵ))*z
    end
    @eval function commute_spin(x::SymbolicUtils.Sym{S1}, y::SymbolicUtils.Sym{S2}) where {S1<:($(Ssym[1])),S2<:($(S[2]))}
        op_x = _to_qumulants(x)
        z = _to_symbolic($(Ssym[3])(op_x.hilbert,op_x.name,op_x.aon))
        y_ = symmetrize(y)
        return x*y_ + y_*x + (im*$(ϵ))*z
    end
    @eval function commute_spin(x::SymbolicUtils.Sym{S1}, y::SymbolicUtils.Sym{S2}) where {S1<:($(S[1])),S2<:($(Ssym[2]))}
        op_x = _to_qumulants(x)
        z = _to_symbolic($(Ssym[3])(op_x.hilbert,op_x.name,op_x.aon))
        x_ = symmetrize(x)
        return x_*y + y*x_ + (im*$(ϵ))*z
    end
end

for (S,Ssym) in zip([:SigmaX,:SigmaY,:SigmaZ],[:SymmetrizedSigmaX,:SymmetrizedSigmaY,:SymmetrizedSigmaZ]) # TODO: cleaner solution for this
    @eval commute_spin(x::SymbolicUtils.Sym{<:($S)},y::SymbolicUtils.Sym{<:($(S))}) = x*y
    @eval commute_spin(x::SymbolicUtils.Sym{<:($S)},y::SymbolicUtils.Sym{<:($(Ssym))}) = x*y
    @eval commute_spin(x::SymbolicUtils.Sym{<:($Ssym)},y::SymbolicUtils.Sym{<:($(S))}) = x*y
end

# for S in [:SigmaX,:SigmaY,:SigmaZ]
#     @eval function rewrite_spinhalf(args_l, args_r, x::SymbolicUtils.Sym{T}, y::SymbolicUtils.Sym{T}) where T<:($(S))
#         if acts_on(x)==acts_on(y)
#             isempty(args_l)&&isempty(args_r)&&(return one(x))
#             return *(args_l...,args_r...)
#         else
#             return nothing
#         end
#     end
# end
#
# for (S,l) in zip(permutations([:SigmaX,:SigmaY,:SigmaZ]), permutations([1,2,3]))
#     ϵ = levicivita(l)
#     @eval function rewrite_spinhalf(args_l, args_r, x::SymbolicUtils.Sym{S1}, y::SymbolicUtils.Sym{S2}) where {S1<:($(S[1])),S2<:($(S[2]))}
#         if acts_on(x)==acts_on(y)
#             op_x = _to_qumulants(x)
#             z = _to_symbolic($(S[3])(op_x.hilbert,op_x.name,op_x.aon))
#             return *(args_l...,($(ϵ)*im)*z,args_r...)
#         else
#             return nothing
#         end
#     end
# end
