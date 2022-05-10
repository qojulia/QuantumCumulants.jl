using LinearAlgebra

"""
Abstract base class for all systems defined in this library.

Currently there are the following concrete systems:

* Spin
* SpinCollection
* CavityMode
* CavitySpinCollection
"""
abstract type System end


"""
A class representing a single spin.

A spin is defined by its position and its detuning to a main
frequency.

# Arguments
* `position`: A vector defining a point in R3.
* `delta`: Detuning.
"""
struct Spin{T1,T2} <: System
    position::T1
    delta::T2
end
Spin(position::Vector; delta::Real=0) = Spin(position, delta)


"""
A class representing a system consisting of many spins.

# Arguments
* `spins`: Vector of spins.
* `polarizations`: Unit vectors defining the directions of the spins.
* `gammas`: Decay rates.
"""
struct SpinCollection{S<:Spin,P<:Vector,G<:Real} <: System
    spins::Vector{S}
    polarizations::Vector{P}
    gammas::Vector{G}
    function SpinCollection{S,P,G}(spins::Vector{S}, polarizations::Vector{P}, gammas::Vector{G}) where {S,P<:Vector{<:Number},G}
        @assert length(polarizations)==length(spins)
        @assert length(gammas)==length(spins)
        new(spins,normalize.(polarizations),gammas)
    end
end
SpinCollection(spins::Vector{S}, polarizations::Vector{P}, gammas::Vector{G}) where {S,P,G} = SpinCollection{S,P,G}(spins, polarizations, gammas)
SpinCollection(spins::Vector{<:Spin}, polarizations::Vector{<:Vector{<:Number}}, gammas::Number) = SpinCollection(spins, polarizations, [gammas for i=1:length(spins)])
SpinCollection(spins::Vector{<:Spin}, polarizations::Vector{<:Number}, args...) = SpinCollection(spins, [polarizations for i=1:length(spins)], args...)
SpinCollection(spins::Vector{<:Spin}, polarizations; gammas=ones(length(spins))) = SpinCollection(spins, polarizations, gammas)

"""
Create a SpinCollection without explicitly creating single spins.

# Arguments
* `positions`: Vector containing the positions of all single spins.
* `polarizations`: Unit vectors defining the directions of the spins.
* `deltas=0`: Detunings.
* `gammas=1`: Decay rates.
"""
function SpinCollection(positions::Vector{<:Vector{<:Real}}, args...; deltas::Union{T,Vector{T}}=zeros(length(positions)), kwargs...) where T<:Real
    if length(deltas)==1
        SpinCollection([Spin(positions[i]; delta=deltas[1]) for i=1:length(positions)], args...; kwargs...)
    else
        SpinCollection([Spin(positions[i]; delta=deltas[i]) for i=1:length(positions)], args...; kwargs...)
    end
end

"""
    interaction.F(ri::Vector, rj::Vector, µi::Vector, µj::Vector)

General F function for arbitrary positions and dipole orientations.

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
"""
function F(ri::Vector, rj::Vector, µi::Vector, µj::Vector)
    rij = ri - rj
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    μi_ = normalize(μi)
    μj_ = normalize(μj)
    T = float(promote_type(eltype(rij),eltype(μi_),eltype(μj_)))
    if rij_norm == 0
        T(2/3.)
    else
        ξ = 2π*rij_norm
        T(dot(µi_, µj_)*(sin(ξ)/ξ + cos(ξ)/ξ^2 - sin(ξ)/ξ^3) + dot(µi_, rijn)*dot(rijn, µj_)*(-sin(ξ)/ξ - 3*cos(ξ)/ξ^2 + 3*sin(ξ)/ξ^3))
    end
end

"""
    interaction.G(ri::Vector, rj::Vector, µi::Vector, µj::Vector)

General G function for arbitrary positions and dipole orientations.

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
"""
function G(ri::Vector, rj::Vector, µi::Vector, µj::Vector)
    rij = ri - rj
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    μi_ = normalize(μi)
    μj_ = normalize(μj)
    T = float(promote_type(eltype(rij),eltype(μi_),eltype(μj_)))
    if rij_norm == 0
        zero(T)
    else
        ξ = 2π*rij_norm
        T(dot(µi_, µj_)*(-cos(ξ)/ξ + sin(ξ)/ξ^2 + cos(ξ)/ξ^3) + dot(µi_, rijn)*dot(rijn, µj_)*(cos(ξ)/ξ - 3*sin(ξ)/ξ^2 - 3*cos(ξ)/ξ^3))
    end
end


"""
    interaction.Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
* γi: Decay rate of first spin.
* γj: Decay rate of second spin.

Note that the dipole moments `μi` and `μj` are normalized internally. To account
for dipole moments with different lengths you need to scale the decay rates
`γi` and `γj`, respectively.
"""
function Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)
    return 0.75*sqrt(γi*γj)*G(ri, rj, µi, µj)
end

"""
    interaction.Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
* γi: Decay rate of first spin.
* γj: Decay rate of second spin.

Note that the dipole moments `μi` and `μj` are normalized internally. To account
for dipole moments with different lengths you need to scale the decay rates
`γi` and `γj`, respectively.
"""
function Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)
    return 1.5*sqrt(γi*γj)*F(ri, rj, µi, µj)
end

"""
    interaction.OmegaMatrix(S::SpinCollection)

Matrix of the dipole-dipole interaction for a given SpinCollection.
"""
function OmegaMatrix(S::SpinCollection)
    spins = S.spins
    mu = S.polarizations
    gamma = S.gammas
    N = length(spins)
    Ω = zeros(Float64, N, N)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        Ω[i,j] = Omega(spins[i].position, spins[j].position, mu[i], mu[j], gamma[i], gamma[j])
    end
    return Ω
end


"""
    interaction.GammaMatrix(S::SpinCollection)

Matrix of the collective decay rate for a given SpinCollection.
"""
function GammaMatrix(S::SpinCollection)
    spins = S.spins
    mu = S.polarizations
    gamma = S.gammas
    N = length(spins)
    return [
        Gamma(spins[i].position, spins[j].position, mu[i], mu[j], gamma[i], gamma[j])
        for i=1:N, j=1:N
    ]
end


"""
    GreenTensor(r::Vector, k::Number=2π)

Calculate the Green's Tensor at position r for wave number k defined by

```math
G = e^{ikr}\\Big[\\left(\\frac{1}{kr} + \\frac{i}{(kr)^2} - \\frac{1}{(kr)^3}\\right)*I -
    \\textbf{r}\\textbf{r}^T\\left(\\frac{1}{kr} + \\frac{3i}{(kr)^2} - \\frac{3}{(kr)^3}\\right)\\Big]
```

Choosing `k=2π` corresponds to the position `r` being given in units of the
wavelength associated with the dipole transition.

Returns a 3×3 complex Matrix.
"""
function GreenTensor(r::Vector{<:Number},k::Real=2π)
    n = norm(r)
    rn = r./n
    return exp(im*k*n)*(
        (1/(k*n) + im/(k*n)^2 - 1/(k*n)^3).*Matrix(I,3,3) +
        -(1/(k*n) + 3im/(k*n)^2 - 3/(k*n)^3).*(rn*rn')
    )
end
