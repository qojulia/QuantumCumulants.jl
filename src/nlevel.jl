"""
    Transition{L,ID,I,J,N,G} <: BasicOperator
    Transition(label::Symbol,i,j,inds;GS=inds[1])
    Transition(label::Symbol,i,j,i_max;GS=inds[1])

Type defining a transition operator from level `j` to `i` in a system consisting
of multiple energy levels.

# Arguments
*`label::Symbol`: The symbolic label of the operator.
*`i`: The index corresponding to the level to which the transition points.
*`j`: The index corresponding from which level the transition occurs.
*`inds`: An iterable object containing all possible indices (levels).

# Optional arguments
*`GS=inds[1]`: The ground state of the system. This specification is used for
    simplification purposes. Since the projectors (transition where `i==j`)
    on each of the levels of the system must form a complete set (sum up to `𝟙`)
    one of them is trivial and can be elminiated. The kwarg `GS` determines
    which operator is eliminated.

# Fields
*`label`: The symbolic label of the operator.
*`id`: An identifier unique to the symbol and indices given to the operator.
    Used for `isequal`.
*`i`: The index corresponding to the level to which the transition points.
*`j`: The index corresponding from which level the transition occurs.
*`inds`: An iterable object containing all possible indices (levels).
*`GS`: The index determining which projector should be eliminated.

# Examples
```jldoctest
julia> σ = Transition(:σ,:g,:e,(:g,:e))
σeg

julia> s = Transition(:s,1,2,2)
s12

julia> s = Transition(:s,1,3,1:3;GS=2)
s13
```
"""
mutable struct Transition{L,ID,I,J,N,G} <: BasicOperator
    label::L
    id::ID
    i::I
    j::J
    inds::N
    GS::G
end

function Transition(label::L,i::I,j::J,inds::N;GS=inds[1]) where {L,I,J,N}
    @assert (i∈inds && j∈inds && GS∈inds)
    return Transition(label,hash.([label,i,j,inds,GS]),i,j,inds,GS)
end
function Transition(label,i::T,j::T,i_max::T;kwargs...) where T<:Int
    return Transition(label,i,j,1:i_max;kwargs...)
end
==(a::Transition,b::Transition) = (a.id==b.id)#(a.label==b.label && a.id==b.id && a.i==b.i && a.j==b.j)
copy(a::Transition) = Transition(a.label,a.id,a.i,a.j,a.inds,a.GS)

Base.adjoint(a::Transition) = Transition(a.label,a.j,a.i,a.inds;GS=a.GS)

function commutation_relation(a::Transition,b::Transition)
    (a.label == b.label && a.inds == b.inds && a.GS==b.GS) || error("Something went wrong here!")
    if a==b
        return zero(a)
    else
        if a.j==b.i
            return Transition(a.label,a.i,b.j,a.inds;GS=a.GS)
        elseif a.i==b.j
            return -Transition(a.label,b.i,a.j,a.inds;GS=a.GS)
        else
            return a*b - b*a
        end
    end
end

function replace_commutator(a::Transition,b::Transition)
    (a.label == b.label && a.inds == b.inds && a.GS==b.GS) || error("Something went wrong here!")
    if a==b
        return (true,simplify_operators(b*a))
    elseif a.j == b.i
        return (true,Transition(a.label,a.i,b.j,a.inds;GS=a.GS))
    else
        return (true,zero(a))
    end
end

ishermitian(a::Transition) = (a.i==a.j)



# Possible solution if we want to be able not to provide indices
# TODO: should we do this?
# import Base: in
# struct EmptyInds end
# in(a,::EmptyInds) = true
# Base.getindex(::EmptyInds, args...) = gensym()
# Base.lastindex(::EmptyInds) = gensym()
# Base.firstindex(::EmptyInds) = gensym()
# Base.isequal(::EmptyInds,::EmptyInds) = true
# Transition(label,i,j) = Transition(label,i,j,EmptyInds())
