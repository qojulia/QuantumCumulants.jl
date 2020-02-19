
"""
    PhotonicOperator <: BasicOperator

Abstract type which is supertype to `Create` and `Destroy`.
"""
abstract type PhotonicOperator <: BasicOperator end

"""
    Destroy <: PhotonicOperator
    Destroy(label::Symbol)

Type defining the photonic annihilation operator.

# Arguments:
*`label::Symbol`: The symbol labelling the operator.

# Fields
*`label`: See above.
*`id`: An identifier unique to the symbol given to the operator. Used for `isequal`.
"""
mutable struct Destroy{L,I} <: PhotonicOperator
    label::L
    id::I
end
Destroy(label::Symbol) = Destroy(label,hash(label))
==(a::Destroy,b::Destroy) = (a.id==b.id)
copy(a::Destroy) = Destroy(a.label,a.id)

"""
    Create <: PhotonicOperator
    Create(label::Symbol)

Type defining the photonic annihilation operator.

# Arguments:
*`label::Symbol`: The symbol labelling the operator.

# Fields
*`label`: See above.
*`id`: An identifier unique to the symbol given to the operator. Used for `isequal`.
"""
mutable struct Create{L,I} <: PhotonicOperator
    label::L
    id::I
end
Create(label::Symbol) = Create(label,hash(label))
==(a::Create,b::Create) = (a.id==b.id)
copy(a::Create) = Create(a.label,a.id)

Base.adjoint(a::Destroy) = Create(a.label,a.id)
Base.adjoint(a::Create) = Destroy(a.label,a.id)

function commutation_relation(a::Destroy,b::Create)
    a.id == b.id || error("Something went wrong here!")
    return one(a)
end
function replace_commutator(a::Destroy,b::Create)
    a.id == b.id || error("Something went wrong here!")
    return (true,b*a + one(a))
end
