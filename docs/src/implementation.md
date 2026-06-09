# Implementation

Let's take a closer look at each step involved from defining a system to arriving at a numerical solution of the underlying time dynamics.


## Hilbert spaces

The first step in treating a system with **QuantumCumulants.jl** is to specify the Hilbert space on which the system is defined. The two most commonly used subspaces are [`FockSpace`](@ref) and [`NLevelSpace`](@ref). The first describes systems whose operators follow the fundamental bosonic commutation relations (such as the quantum harmonic oscillator), whereas the latter describes systems consisting of a finite number of energy levels with an arbitrary energy difference in between (such as atoms). Additional subspaces ([`PauliSpace`](@ref), [`SpinSpace`](@ref), [`PhaseSpace`](@ref)) cover spin systems and position/momentum.

A [`FockSpace`](@ref) simply needs a name in order to be defined:

```@example hilbert-space
using QuantumCumulants # hide
hf = FockSpace(:fock1)
nothing # hide
```

[`NLevelSpace`](@ref) requires a name as well as labels for the energy levels. For example

```@example hilbert-space
h_atom = NLevelSpace(:atom, (:g, :e))
nothing # hide
```

defines an [`NLevelSpace`](@ref) with the name `:atom` and the two levels labeled by `:g` and `:e`, respectively. The levels can be labeled by (almost) anything. For example, `NLevelSpace(:two_level, (1,2))` would define a Hilbert space describing a system with the two discrete energy levels labeled by `1` and `2`. Specifically for numbers, there is also the short-hand method to write `NLevelSpace(:five_level, 5)` which creates a system with levels `1:5`. By default the first level in the list of all levels is designated as the ground state. This can be changed by specifying the ground state explicitly as a third argument to [`NLevelSpace`](@ref), e.g. `NLevelSpace(:four_level, 4, 2)` would designate the state `2` as the ground state. The ground state projector will be eliminated during simplification (see below).

Composite systems are generally described by a [`ProductSpace`](@ref), i.e. a Hilbert space that consists of multiple subspaces. They can be created using the [`tensor`](@ref) function or the unicode symbol [`⊗`](@ref) [\otimes]. For example

```@example hilbert-space
h_prod1 = tensor(hf, h_atom)
h_prod2 = tensor(h_prod1, NLevelSpace(:three_level, 3))
h_prod3 = tensor(hf, h_atom, NLevelSpace(:three_level, 3)) # == h_prod2
nothing # hide
```

creates two product spaces. The first, `h_prod1`, consists of the previously defined `FockSpace(:fock1)` and `NLevelSpace(:atom, (:g,:e))`. The second one, `h_prod2`, adds in another `NLevelSpace(:three_level, 3)`. In principle arbitrarily many systems can be combined this way.


## Operators (a.k.a. *q*-numbers)

Once the Hilbert space of the system has been defined, we can proceed by defining operators, or *q*-numbers, on them. They are the fundamental building blocks of symbolic expressions in **QuantumCumulants.jl**. The two most common operators are the quantum harmonic destruction operator [`Destroy`](@ref) which acts on a [`FockSpace`](@ref), and a [`Transition`](@ref) operator which describes a transition between any two levels on an [`NLevelSpace`](@ref). These operators can only be defined on the corresponding Hilbert spaces.
Other built-in *q*-number types ([`Pauli`](@ref), [`Spin`](@ref), [`Position`](@ref), [`Momentum`](@ref)) are also available, and custom operator types can be added (see the [interface](@ref interface) section below).

Here are a few examples:

```@example operators
using QuantumCumulants # hide
hf = FockSpace(:fock)
a = Destroy(hf, :a)

h_atom = NLevelSpace(:atom, (:g, :e))
σge = Transition(h_atom, :σ, :g, :e)
σ(i, j) = Transition(h_atom, :σ, i, j)
@assert isequal(σge, σ(:g, :e)) # true
nothing # hide
```

As you can see, the destruction operator [`Destroy`](@ref) is created on a [`FockSpace`](@ref) and given a name. The transition operator additionally requires you to specify the levels between which it describes the transition. Wrapping the constructor in a small `σ(i, j) = Transition(...)` closure gives the familiar callable shape used throughout the docs. In bra-ket notation, the transition operator `Transition(h, i, j)` is simply ``|i\rangle \langle j|``. The bosonic creation operator is simply given by the `adjoint` of [`Destroy`](@ref).

These fundamental operators are all subtypes of [`QSym`](@ref), and constitute the basic symbolic building blocks for the noncommutative algebra used in **QuantumCumulants.jl**. They can be combined using standard algebraic functions.

```@example operators
ex_fock = 0.1*a'*a
ex_trans = im*(σ(:g, :e) - σ(:e, :g))
nothing # hide
```

Note that only operators that are defined on the same Hilbert space can be algebraically combined. The resulting expressions are stored as [`QAdd`](@ref) (a sum of products with prefactors). Internally there is no separate `QMul` value: products are stored as `QTerm` dict-keys inside a `QAdd`, so every non-atomic expression normalises to a `QAdd`.

In composite systems, we also need to specify on which subsystem the respective operator acts. This information is important as operators acting on different subsystems commute with one another, but operators acting on the same one do not. When multiplying together operators in a composite system, they are automatically ordered according to the order of Hilbert spaces. It's specified by an additional argument when creating operators.

```@example operators
h_prod = FockSpace(:fock1) ⊗ FockSpace(:fock2)
a = Destroy(h_prod, :a, 1)
b = Destroy(h_prod, :b, 2)
a*b # a*b
b*a # a*b
a'*b*a # a'*a*b
nothing # hide
```

If a subspace occurs only once in a [`ProductSpace`](@ref), the choice on which an operator acts is unique and can therefore be omitted on construction.

```@example operators
h_prod = FockSpace(:fock1) ⊗ FockSpace(:fock2) ⊗ NLevelSpace(:atom, (:g, :e))
σ(i, j) = Transition(h_prod, :σ, i, j) # no need to specify acts_on
nothing # hide
```

For convenience, there is also a macro that can be used to construct operators:

```@example operators
h = FockSpace(:fock) ⊗ NLevelSpace(:two_level, 2)
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)
ex = a'*σ(1, 2) + a*σ(2, 1)
nothing # hide
```

## Symbolic parameters

Symbolic scalar parameters are declared with `@variables` from [**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl), which is re-exported via [**SecondQuantizedAlgebra.jl**](https://github.com/qojulia/SecondQuantizedAlgebra.jl). The default symbolic number type is complex; constrain to a real number by adding `::Real`. Symbolic parameters compose with *q*-numbers to build Hamiltonians, e.g.

```@example c-numbers
using QuantumCumulants # hide
h = FockSpace(:fock)
@variables ω::Real η::Real
@qnumbers a::Destroy(h)
H = ω*a'*a + η*(a + a')
nothing # hide
```

For indexed (per-site) parameters use [`IndexedVariable`](@ref) and [`DoubleIndexedVariable`](@ref) (see [Symbolic Sums and Indices](@ref)).


## Operator expressions and commutation relations

The equations of motion of *q*-numbers are determined by evaluating commutators. This can be done by using fundamental commutation relations, which are immediately applied whenever operators are combined in an algebraic expression.

For the quantum harmonic oscillator destruction operator ``a``, we have the canonical commutator

```math
[a,a^\dagger] = 1.
```

Within the framework, we choose normal ordering, which surmounts to the rewriting rule

```math
a a^\dagger ~\Rightarrow~ a^\dagger a +1.
```

For transition operators ``\sigma^{ij}`` denoting a transition from level ``j`` to level ``i``, on the other hand, we have a rule for products,

```math
\sigma^{ij}\sigma^{kl} ~\Rightarrow~ \delta_{jk}\sigma^{il},
```

which is implemented as a rewriting rule just so. Additionally, we use the fact that in a system with levels ``\{1,...,n\}``

```math
\sum_{j=1}^n \sigma^{jj} = 1
```

in order to eliminate the projector on the ground state. This reduces the amount of equations required for each [`NLevelSpace`](@ref) by 1. As mentioned before, the ground state is by default chosen to be the first (but this can be changed). Hence, the default rewriting rule to eliminate the ground-state projector is

```math
\sigma^{11} ~\Rightarrow~ 1 - \sum_{j=2}^n \sigma^{jj}.
```

Any expression involving operators is stored as a [`QAdd`](@ref). The expression tree is structured so that the application of commutation relations can be done efficiently: products are stored as keys in a `QTermDict` mapping each operator product to its prefactor. This makes it easy and efficient to recurse through the expression and find pairs of operators that should be rewritten according to some commutation relation.

Note that only simplification using commutation relations is implemented directly in **QuantumCumulants.jl**. For any other simplification routines, operators are averaged (without applying a cumulant expansion) which makes them numbers. Those numbers are stored as expressions in [**SymbolicUtils.jl**](https://github.com/JuliaSymbolics/SymbolicUtils.jl) and simplified according to standard simplification rules. Afterwards, they can be converted back into operator expressions with `undo_average`.

Here's a short example:

```@example meanfield
using QuantumCumulants # hide
h = FockSpace(:fock)
@qnumbers a::Destroy(h)
a*a' # returns a'*a + 1
nothing # hide
```

In order to derive equations of motion, you need to specify a Hamiltonian and the operator (or a list of operators) of which you want to derive the Heisenberg equations and pass them to [`meanfield`](@ref), which stores both the operator equations as well as the averaged ones. In the end, we only want to work with averages.

```@example meanfield
using Latexify # hide
set_default(double_linebreak=true) # hide
@variables ω::Real η::Real
H = ω*a'*a + η*(a + a') # Driven cavity Hamiltonian
me = meanfield([a, a'*a], H)
```

## Cumulant expansion

Averaging (using [`average`](@ref)) and the [`cumulant_expansion`](@ref) are essential to convert the system of *q*-number equations to *c*-number equations. Averaging alone converts any operator product to a *c*-number, yet you will not arrive at a closed set of equations without truncating at a specific order. An average is stored as a symbolic expression: each `Average` carries the underlying operator in its metadata so the original `QAdd` can be recovered with `undo_average`.

The order of an average is given by the number of constituents in the product. For example

```@example cumulant
using QuantumCumulants # hide
h = FockSpace(:fock)
@qnumbers a::Destroy(h)

avg1 = average(a)
get_order(avg1) # 1

avg2 = average(a'*a)
get_order(avg2) # 2
nothing # hide
```

The cumulant expansion can then be used to express any average by averages up to a specified order (see also the [theory section](@ref theory)):

```@example cumulant
cumulant_expansion(avg2, 1)
nothing #hide
```

When deriving the equations of motion using the [`meanfield`](@ref) function, the [`cumulant_expansion`](@ref) is immediately applied if you specify an order, e.g. `meanfield(ops, H; order=2)`. Before you can actually solve the system of equations, you need to ensure that it is complete, i.e. there are no averages missing. This can be checked with [`find_missing`](@ref). Alternatively, you can automatically complete a system of equations using the [`complete`](@ref) function which will internally use [`find_missing`](@ref) to look for missing averages and derive equations for those.


## Numerical solution

To solve a system of equations, we convert a set of [`MeanFieldEquations`](@ref) to a [`ModelingToolkitBase.System`](https://github.com/SciML/ModelingToolkitBase.jl), which represents a symbolic set of ordinary differential equations. The resulting system can be compiled into fast numerical functions that feed directly into the [**OrdinaryDiffEq.jl**](https://github.com/SciML/OrdinaryDiffEq.jl) solvers. ModelingToolkitBase also offers symbolic Jacobians for better stability and performance.

To obtain a `System` from a [`MeanFieldEquations`](@ref), call `System` and then `mtkcompile`:

```@example meanfield
using ModelingToolkitBase
sys = mtkcompile(System(me; name=:cavity))
nothing # hide
```

To obtain a numerical solution we construct an `ODEProblem` and solve it.

```@example meanfield
using OrdinaryDiffEq
u0 = zeros(ComplexF64, length(me.states))
p0 = [ω => 1.0, η => 0.1]
prob = ODEProblem(sys, merge(initial_values(me, u0), Dict(p0)), (0.0, 1.0))
sol = solve(prob, Tsit5())
nothing # hide
```

The state of the system at each time-step is stored in `sol.u`. To extract the trajectory of a specific average we use [`get_solution`](@ref), which substitutes the symbolic average into the compiled solution and returns a callable in `t`. This works both for averages that sit on the LHS of `me` and for derived products that don't:

```@example meanfield
ts  = range(0.0, 1.0; length=100)
a_t = get_solution(sol, a, me).(ts)
n_t = real.(get_solution(sol, a'*a, me).(ts))
nothing # hide
```


### Calculating the initial state

When trying to solve a system of equations numerically, it can sometimes become tricky to find the correct initial state.
In the above, we simply did `u0 = zeros(ComplexF64, length(me.states))`, since that was a viable initial state.
However, things become more involved when you have, say, a superposition of two coherent states in a harmonic oscillator as starting point,

```math
|\psi_0\rangle = \frac{1}{\sqrt{2}} \left( |\alpha\rangle + |\beta\rangle \right),
```

where $\alpha, \beta \in \mathbb{C}$ are the respective complex amplitudes.
While computing the first-order expectation values such as `\langle \psi_0 | a |\psi_0\rangle` is still simple enough, things become more tricky in higher orders and when mixing in another Hilbert space (e.g. an atom in a cavity).
Since the system of equations can become quite large, this may result in quite some manual effort when trying to calculate all initial values.
And we hate manual effort.

These expectations values are, however, only difficult to calculate symbolically, yet are easy enough to compute numerically.
**QuantumCumulants.jl** therefore offers a convenient integration to [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl), which allows you to quickly calculate the initial expectation values of a system of equations from a given numerical initial state.
The function is called [`initial_values`](@ref).
For example, we could use it in the above to compute a coherent initial state:

```@example meanfield
using QuantumOpticsBase
b = FockBasis(10)
alpha = 0.3 + 0.4im
psi_0 = coherentstate(b, alpha)
u0 = initial_values(me, psi_0)
nothing # hide
```

You can also compute initial values for mixed states by passing a density operator:

```@example meanfield
u0 = initial_values(me, dm(psi_0))
nothing # hide
```

#### Levels for `NLevelSpace`

The conversion to a numeric representation between [`FockSpace`](@ref) and `FockBasis` is always uniquely defined.
For an [`NLevelSpace`](@ref) the mapping from symbolic levels to the basis states of `NLevelBasis` is fixed by the level labels you choose at construction time: integer labels map directly to the matching basis index.

```@example levelmap
using QuantumCumulants, QuantumOpticsBase
h = NLevelSpace(:TwoLevelAtom, 2)   # integer levels 1:2
b = NLevelBasis(2)
s = Transition(h, :s, 1, 2)
@assert to_numeric(s, b) == transition(b, 1, 2)
nothing # hide
```

To numerically convert an [`NLevelSpace`](@ref) defined with symbolic labels (e.g. `(:g, :e)`), declare the space with integer labels in the first place, since the `level_map` keyword from the 0.4 series is no longer accepted.

#### Numeric averages and conversion

While the examples so far were relatively simple and would have been easy to calculate by hand, things quickly become more difficult whenever product spaces and higher-order products are involved.

Behind the scenes, [`initial_values`](@ref) just uses the [`numeric_average`](@ref) method in order to compute the numeric expectation value for the given operators and states.
This method in turn calls into the numeric conversion [`to_numeric`](@ref) and then uses `QuantumOpticsBase.expect` on the result in order to calculate the respective expectation values for the given state and operators numerically.
Should you need to compute numerical averages from a symbolic one for a given numerical state, you can also call [`numeric_average`](@ref) directly.

```@example tonumeric
using QuantumCumulants, QuantumOpticsBase
hfock   = FockSpace(:cavity)
hnlevel = NLevelSpace(:ThreeLevelAtom, 3)   # integer levels 1:3
h       = hfock ⊗ hnlevel
a       = Destroy(h, :a)
s       = Transition(h, :s, 3, 1)           # |3⟩⟨1|

bfock   = FockBasis(10)
bnlevel = NLevelBasis(3)
psi     = coherentstate(bfock, 0.3) ⊗
          (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) / sqrt(2)

avg     = average(a' * s)
avg_num = numeric_average(avg, psi)
nothing # hide
```

Similarly, you can also just obtain the numerical representation of an operator by directly calling [`to_numeric`](@ref) on a given basis.

```@example tonumeric
b = bfock ⊗ bnlevel
a_num = to_numeric(a, b)
nothing # hide
```

Note that [`to_numeric`](@ref) returns a `SparseOperator` for single operators, but a `LazyTensor` operator whenever a product space is involved.
Lazy evaluation of tensor products is incredibly useful here, as symbolically easy to treat systems can become quite large numerically.

When a large number of Hilbert spaces is involved, it can even become tricky to store a single `Ket`.
In order to overcome this limitation, QuantumOpticsBase also offers lazy evaluation of state products, allowing you to compute expectation values and initial states for very large product states.

```@example tonumeric
psi_lazy = LazyKet(b, (
    coherentstate(bfock, 0.3),
    (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) / sqrt(2),
))
avg_num_lazy = numeric_average(avg, psi_lazy)
@assert isapprox(avg_num, avg_num_lazy)
```


## [The *q*-number interface](@id interface)

The built-in operator types ([`Destroy`](@ref), [`Transition`](@ref), [`Pauli`](@ref), [`Spin`](@ref), [`Position`](@ref), [`Momentum`](@ref)) all derive from `SecondQuantizedAlgebra.QSym` and talk to the algebra through a small set of hooks. Adding a new operator type requires:

* a struct subtyping `SecondQuantizedAlgebra.QSym` whose fields encode the operator's identity. The convention used by the built-ins is `name::Symbol`, `space_index::Int`, `index::SecondQuantizedAlgebra.Index` (use `SecondQuantizedAlgebra.NO_INDEX` when there's no symbolic index), plus any operator-specific fields,
* `Base.isequal`, `Base.:(==)`, and `Base.hash` over all identity-carrying fields,
* `Base.adjoint`, even if the operator is Hermitian (return `op`),
* the canonicalisation hooks consumed by SQA's streaming sorter: `_site_compare`, `_can_commute`, `_commute_pair`, `_reduce_pair` (and `_ground_state_expand` for level-projector logic). SQA provides trivial cross-type fallbacks, so you only need methods between the new type and itself (and any other type it shares a site with),
* a `_type_order` slot (both `_type_order(::Type{MyOp})` and `_type_order(::MyOp)`) so the canonical sort can break ties between different operator types.

#### Example: Harmonic oscillator quadratures

To illustrate, we re-implement the harmonic-oscillator position ``x`` and momentum ``p`` as a custom pair. (**QuantumCumulants.jl** already ships [`Position`](@ref) / [`Momentum`](@ref) on [`PhaseSpace`](@ref); the point of this example is to walk through the hooks.) The pair satisfies

```math
[p, x] = -i \quad \Rightarrow \quad p\,x = x\,p - i.
```

The two operators live on a `MyPhaseSpace`:

```@example custom-operators
using Latexify # hide
set_default(double_linebreak=true) # hide
using QuantumCumulants
using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

struct MyPhaseSpace <: SQA.HilbertSpace
    name::Symbol
end
Base.:(==)(a::MyPhaseSpace, b::MyPhaseSpace) = a.name == b.name
Base.hash(a::MyPhaseSpace, h::UInt) = hash(:MyPhaseSpace, hash(a.name, h))

struct MyPosition <: SQA.QSym
    name::Symbol
    space_index::Int
    index::SQA.Index
end
MyPosition(name::Symbol, si::Int) = MyPosition(name, si, SQA.NO_INDEX)
MyPosition(h::MyPhaseSpace, name::Symbol) = MyPosition(name, 1)

struct MyMomentum <: SQA.QSym
    name::Symbol
    space_index::Int
    index::SQA.Index
end
MyMomentum(name::Symbol, si::Int) = MyMomentum(name, si, SQA.NO_INDEX)
MyMomentum(h::MyPhaseSpace, name::Symbol) = MyMomentum(name, 1)
nothing # hide
```

Equality, hashing, and adjoint are mechanical:

```@example custom-operators
Base.adjoint(op::MyPosition) = op
Base.adjoint(op::MyMomentum) = op

Base.isequal(a::MyPosition, b::MyPosition) =
    a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.isequal(a::MyMomentum, b::MyMomentum) =
    a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::MyPosition, b::MyPosition) = isequal(a, b)
Base.:(==)(a::MyMomentum, b::MyMomentum) = isequal(a, b)

Base.hash(a::MyPosition, h::UInt) =
    hash(:MyPosition, hash(a.name, hash(a.space_index, hash(a.index, h))))
Base.hash(a::MyMomentum, h::UInt) =
    hash(:MyMomentum, hash(a.name, hash(a.space_index, hash(a.index, h))))
nothing # hide
```

The five operator hooks. `_site_compare` determines a partial order over operators on the same site (for the canonical sort), `_can_commute` says whether two same-site operators swap with no residual, `_commute_pair` returns `(swapped_b, swapped_a, residual_coefficient, residual_ops)` when they don't, and `_reduce_pair` collapses local algebraic identities (we have none here, so we let the fallback handle it):

```@example custom-operators
function SQA._site_compare(a::MyPosition, b::MyPosition, ne::Vector{SQA.NonEqualPair})
    a.space_index == b.space_index ||
        return a.space_index < b.space_index ? SQA.Less : SQA.Greater
    a.name == b.name || return a.name < b.name ? SQA.Less : SQA.Greater
    return a.index == b.index ? SQA.Equal : SQA.Undetermined
end
function SQA._site_compare(a::MyMomentum, b::MyMomentum, ne::Vector{SQA.NonEqualPair})
    return SQA._site_compare(
        MyPosition(a.name, a.space_index, a.index),
        MyPosition(b.name, b.space_index, b.index), ne,
    )
end
function SQA._site_compare(a::MyPosition, b::MyMomentum, ne::Vector{SQA.NonEqualPair})
    a.space_index == b.space_index ||
        return a.space_index < b.space_index ? SQA.Less : SQA.Greater
    return a.index == b.index ? SQA.Equal : SQA.Undetermined
end
function SQA._site_compare(a::MyMomentum, b::MyPosition, ne::Vector{SQA.NonEqualPair})
    a.space_index == b.space_index ||
        return a.space_index < b.space_index ? SQA.Less : SQA.Greater
    return a.index == b.index ? SQA.Equal : SQA.Undetermined
end

# Canonical order: x before p on the same site. x·x and p·p commute trivially.
SQA._can_commute(::MyPosition, ::MyPosition) = true
SQA._can_commute(::MyMomentum, ::MyMomentum) = true
SQA._can_commute(::MyPosition, ::MyMomentum) = true   # already in canonical order
SQA._can_commute(::MyMomentum, ::MyPosition) = false  # needs the residual swap

# p · x = x · p - i · I
SQA._commute_pair(p::MyMomentum, x::MyPosition) = (x, p, SQA._to_cnum(-im), SQA._EMPTY_OPS)

# Order slot for the canonical sort. Pick any integers not used by SQA's
# built-ins (Destroy=0 … Momentum=6); the absolute values don't matter, only
# the relative order between the new types and any others they coexist with.
SQA._type_order(::Type{MyPosition}) = 100
SQA._type_order(::Type{MyMomentum}) = 101
SQA._type_order(::MyPosition) = 100
SQA._type_order(::MyMomentum) = 101
nothing # hide
```

Optionally, teach Latexify how to render the new operators so equations display nicely. SQA ships a `@latexrecipe` per operator type (`Destroy`, `Transition`, `Position`, …) with no generic fallback, so a custom type needs its own:

```@example custom-operators
@latexrecipe function f(x::MyPosition)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end
@latexrecipe function f(x::MyMomentum)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end
nothing # hide
```

That's enough machinery for the canonicalisation pipeline. We can now build a Hamiltonian and derive equations of motion:

```@example custom-operators
h = MyPhaseSpace(:osc)
x = MyPosition(h, :x)
p = MyMomentum(h, :p)

@variables ω::Real m::Real
H = p^2/(2m) + 0.5*m*ω^2*x^2

eqs = meanfield([x, p], H)
```

The pattern generalises to any new operator family: copy the field layout (`name`, `space_index`, `index`, plus operator-specific fields), define the five hooks against the existing operator types you need to interact with, and let SQA's pipelines handle products, sorting, and indexed sums automatically. The built-in `Pauli`, `Spin`, `Position`, `Momentum`, and `Transition` definitions under `SecondQuantizedAlgebra.jl/src/operators/` are the closest reference implementations.
