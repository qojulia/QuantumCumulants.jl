# Implementation

Let's take a closer look at each step involved from defining a system to arriving at a numerical solution of the underlying time dynamics.


## Hilbert spaces

The first step in treating a system with **QuantumCumulants.jl** is to specify the Hilbert space on which the system is defined. There are two types of Hilbert spaces implemented, namely [`FockSpace`](@ref) and [`NLevelSpace`](@ref). The first describes systems whose operators follow the fundamental bosonic commutation relations (such as the quantum harmonic oscillator), whereas the latter describes systems consisting of a finite number of energy levels with an arbitrary energy difference in between (such as atoms).

A [`FockSpace`](@ref) simply needs a name in order to be defined:

```@example hilbert-space
using QuantumCumulants # hide
hf = FockSpace(:fock1)
nothing # hide
```

[`NLevelSpace`](@ref) requires a name as well as labels for the energy levels. For example

```@example hilbert-space
h_atom = NLevelSpace(:atom, (:g,:e))
nothing # hide
```

defines an [`NLevelSpace`](@ref) with the name `:atom` and the two levels labeled by `:g` and `:e`, respectively. Note that the levels can be labeled by (almost) anything. For example, `NLevelSpace(:two_level, (1,2))` would define a Hilbert space describing a system with the two discrete energy levels labeled by `1` and `2`. Specifically for numbers, there is also the short-hand method to write `NLevelSpace(:five_level, 5)` which creates a system with levels `1:5`. Note that by default the first level in the list of all levels is designated as the ground state. This can be changed by specifying the ground state explicitly as a third argument to [`NLevelSpace`](@ref), e.g. `NLevelSpace(:four_level, 4, 2)` would designate the state `2` as the ground state. The ground state projector will be eliminated during simplification (see below).

Composite systems are generally described by a [`ProductSpace`](@ref), i.e. a Hilbert space that consists of multiple subspaces. Each subspace is either a [`FockSpace`](@ref) or an [`NLevelSpace`](@ref). They can be created using the [`tensor`](@ref) function or the unicode symbol [`⊗`](@ref) [\otimes]. For example

```@example hilbert-space
h_prod1 = tensor(hf, h_atom)
h_prod2 = tensor(h_prod1, NLevelSpace(:three_level, 3))
h_prod3 = tensor(hf, h_atom, NLevelSpace(:three_level, 3)) # == h_prod2
nothing # hide
```

creates two product spaces. The first, `h_prod1`, consists of the previously defined `FockSpace(:fock1)` and `NLevelSpace(:atom, (:g,:e))`. The second one, `h_prod2`, adds in another `NLevelSpace(:three_level, 3)`. In principle arbitrarily many systems can be combined this way.


## Operators (a.k.a. *q*-numbers)

Once the Hilbert space of the system has been defined, we can proceed by defining operators, or *q*-numbers, on them. They are the fundamental building blocks of symbolic expressions in **QuantumCumulants.jl**. Again, there are essentially two kinds of operators implemented: the quantum harmonic destruction operator [`Destroy`](@ref) which acts on a [`FockSpace`](@ref), as well as a [`Transition`](@ref) operator which describes a transition between any two levels on an [`NLevelSpace`](@ref). These operators can only be defined on the corresponding Hilbert spaces.
Note that there is no intrinsic reason that prevents us from implementing more types of operators ([see below](@ref interface)), there was simply no need to do that so far.

Here are a few examples:

```@example operators
using QuantumCumulants # hide
hf = FockSpace(:fock)
a = Destroy(hf, :a)

h_atom = NLevelSpace(:atom,(:g,:e))
σge = Transition(h_atom, :σ, :g, :e)
σ = Transition(h_atom, :σ)
@assert isequal(σge, σ(:g,:e)) # true
nothing # hide
```

As you can see, the destruction operator [`Destroy`](@ref) is created on a [`FockSpace`](@ref) and given a name. The transition operator, however, additionally requires you to specify the levels between which it describes the transition. Defining a transition without levels specified creates a callable instance which needs to be called with valid level labels before one can actually use it in any algebraic expressions. Note that in Bra-Ket notation, the transition operator `Transition(h, i, j)` is simply ``|i\rangle \langle j|``. Also, the bosonic creation operator is simply given by the `adjoint` of [`Destroy`](@ref).

These fundamental operators are all of subtypes of [`QSym`](@ref), and constitute the basic symbolic building blocks for the noncommutative algebra used in **QuantumCumulants.jl**. They can be combined using standard algebraic functions.

```@example operators
ex_fock = 0.1*a'*a
ex_trans = im*(σ(:g,:e) - σ(:e,:g))
nothing # hide
```

Note that only operators that are defined on the same Hilbert space can be algebraically combined. The resulting expressions are stored as [`QTerm`](@ref) types.

In composite systems, we also need to specify on which subsystem the respective operator acts. This information is important as operators acting on different subsystems commute with one another, but operators acting on the same one do not. When multiplying together operators in a composite systems, they are automatically ordered according to the order of Hilbert spaces. It's specified by an additional argument when creating operators.

```@example operators
h_prod = FockSpace(:fock1) ⊗ FockSpace(:fock2)
a = Destroy(h_prod,:a,1)
b = Destroy(h_prod,:b,2)
a*b # a*b
b*a # a*b
a'*b*a # a'*a*b
nothing # hide
```

If a subspace occurs only once in a [`ProductSpace`](@ref), the choice on which an operator acts is unique and can therefore be omitted on construction.

```@example operators
h_prod = FockSpace(:fock1) ⊗ FockSpace(:fock2) ⊗ NLevelSpace(:atom,(:g,:e))
σ = Transition(h_prod, :σ) # no need to specify acts_on
nothing # hide
```

For convenience, there is also a macro that can be used to construct operators:

```@example operators
h = FockSpace(:fock) ⊗ NLevelSpace(:two_level, 2)
@qnumbers a::Destroy(h) σ::Transition(h)
ex = a'*σ(1,2) + a*σ(2,1)
nothing # hide
```

## Symbolic parameters (a.k.a. *c*-numbers)

Commutative numbers (*c*-numbers) are represented by `SymbolicUtils.Sym` from the [**SymbolicUtils.jl**](https://github.com/JuliaSymbolics/SymbolicUtils.jl) package and a custom subtype to `Number` called [`CNumber`](@ref). They are generally assumed to be complex numbers and can be defined with the [`cnumbers`](@ref) function or the corresponding macro [`@cnumbers`](@ref). You can use them together with *q*-numbers to build symbolic expressions describing the Hamiltonian, e.g.

```@example c-numbers
using QuantumCumulants # hide
h = FockSpace(:fock)
@cnumbers ω η
@qnumbers a::Destroy(h)
H = ω*a'*a + η*(a + a')
nothing # hide
```

Real numbers (*r*-numbers) are similar to *c*-numbers, except that they are their own complex conjugate. They can be defined with the [`rnumbers`](@ref) function or the corresponding macro [`@rnumbers`](@ref). 

```@example r-numbers
using QuantumCumulants # hide
@rnumbers ω η
ω' # ω 
exp(1im*η)*(exp(1im*η))' # 1
nothing # hide
```


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

which is implemented as rewriting rule just so. Additionally, we use the fact that in a system with levels ``\{1,...,n\}``

```math
\sum_{j=1}^n \sigma^{jj} = 1
```

in order to eliminate the projector on the ground state. This reduces the amount of equations required for each [`NLevelSpace`](@ref) by 1. Note that, as mentioned before, the ground state is by default chosen to be the first (but this can be changed). Hence, the default rewriting rule to eliminate the ground-state projector is

```math
\sigma^{11} ~\Rightarrow~ 1 - \sum_{j=2}^n \sigma^{jj}.
```

Any expression involving operators is stored as a [`QTerm`](@ref) type. The expression trees are structured such that the application of commutation relations can be done efficiently. There are two concrete subtypes of [`QTerm`](@ref), namely `QMul` representing a multiplication and `QAdd` representing an addition. Methods of multiplication and addition are implemented such that `QSym < QMul < QAdd`, i.e. a multiplication can only consist of numbers and fundamental operators (it cannot contain another multiplication or addition) and `QAdd` is always at the highest level possibly containing numbers, `QSym`s and `QMul`s (but no other `QAdd`s). This makes it easy and efficient to recurse through the expression tree and find pairs of operators that should be rewritten according to some commutation relation.

Note that only simplification using commutation relations is implemented directly in **QuantumCumulants.jl**. For any other simplification routines, operators are averaged (without applying a cumulant expansion) which makes them numbers. Those numbers are stored as expressions in [**SymbolicUtils.jl**](https://github.com/JuliaSymbolics/SymbolicUtils.jl) and simplified according to standard simplification rules. Afterwards, they can be converted back into [`QTerm`](@ref) expressions.

Here's a short example:

```@example meanfield
using QuantumCumulants # hide
h = FockSpace(:fock)
@qnumbers a::Destroy(h)
a*a' # returns a'*a + 1
nothing # hide
```

In order to derive equations of motion, you need to specify a Hamiltonian and the operator (or a list of operators) of which you want to derive the Heisenberg equations and pass them to [`meanfield`](@ref), which stores both the operator as well as the average equations. In the end, we only want to work with averages.

```@example meanfield
using Latexify # hide
set_default(double_linebreak=true) # hide
@cnumbers ω η
H = ω*a'*a + η*(a + a') # Driven cavity Hamiltonian
me = meanfield([a, a'*a], H)
```

## Cumulant expansion

Averaging (using [`average`](@ref)) and the [`cumulant_expansion`](@ref) are essential to convert the system of *q*-number equations to *c*-number equations. Averaging alone converts any operator product to a *c*-number, yet you will not arrive at a closed set of equations without truncating at a specific order. An average is stored as a symbolic expression. Specifically, the average of an operator `op` is internally represented by `SymbolicUtils.Term{AvgSym}(sym_average, [op])`.

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
average(a'*a, 1) # short-hand for cumulant_expansion(average(a'*a), 1)
nothing #hide
```

When deriving the equations of motion using the [`meanfield`](@ref) function, the [`cumulant_expansion`](@ref) is immediately applied if you specify an order, e.g `meafield(ops,H;order=2)`.
Before you can actually solve the system of equations, you need to ensure that it is complete, i.e. there are no averages missing. This can be checked with [`find_missing`](@ref). Alternatively, you can automatically complete a system of equations using the [`complete`](@ref) function which will internally use [`find_missing`](@ref) to look for missing averages and derive equations for those.


## Numerical solution

Finally, in order to actually solve a system of equations, we need to convert a set of equations to an [`System`](https://docs.sciml.ai/ModelingToolkit/dev/API/System/), which represents a symbolic set of ordinary differential equations. `System`s are part of the [**ModelingToolkit.jl**](https://github.com/SciML/ModelingToolkit.jl) framework, which allows for generating fast numerical functions that can be directly used in the [**OrdinaryDiffEq.jl**](https://github.com/SciML/OrdinaryDiffEq.jl) package. On top of that, [**ModelingToolkit.jl**](https://github.com/SciML/ModelingToolkit.jl) also offers a variety of additional functionality, such as the symbolic computation of Jacobians for better stability and performance.

To obtain an `System` from [`MeanfieldEquations`](@ref), you simply need to call the constructor:

```@example meanfield
using ModelingToolkit
@named sys = System(me)
nothing # hide
```

Finally, to obtain a numerical solution we can construct an `ODEProblem` and solve it.

```@example meanfield
using OrdinaryDiffEq
p0 = Dict(ω => 1.0, η => 0.1)
u0 = Dict(zip(unknowns(sys),zeros(ComplexF64, length(me))))
prob = ODEProblem(sys,merge(u0,p0),(0.0,1.0))
sol = solve(prob, RK4())
nothing # hide
```

Now, the state of the system at each time-step is stored in `sol.u`. To access one specific solution, you can simply type e.g. `sol[average(a)]` to obtain the time evolution of the expectation value ``\langle a \rangle``.


### Calculating the initial state

When trying to solve a system of equations numerically, it can sometimes become tricky to find the correct initial state.
In the above, we simply did `u0 = zeros(ComplexF64, length(me))`, since that was a viable initial state.
However, things become more involved when you have, say, a superposition of two coherent states in a harmonic oscillator as starting point,

```math
|\psi_0\rangle = \frac{1}{\sqrt{2}} \left( |\alpha\rangle + |\beta\rangle \right),
```

where $\alpha, \beta \in \mathbb{C}$ are the respective complex amplitudes.
While computing the first-order expectation values such as `\langle \psi_0 | a |\psi_0\rangle` is still simple enough, things become more tricky in higher orders and when mixing in another Hilbert space (e.g. an atom in a cavity).
Since the system of equations can become quite large, this may result in quite some manual effort when trying to calculate all initial values.
And we hate manual effort.

These expectations values are, however, only difficult to calculate symbolically, yet are easy enough to compute numerically.
QuantumCumulants therefore offers a convenient integration to [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl), which allows you to quickly calculate the initial expectation values of a system of equations from a given numerical initial state.
The function is called [`initial_values`](@ref).
For example, we could use it in the above example to compute a coherent initial states

```@example meanfield
using QuantumOpticsBase
b = FockBasis(10)
alpha = 0.3 + 0.4im
psi_0 = coherentstate(b, alpha)
u0 = initial_values(me, psi_0)
nothing # hide
```

Note that you can also compute initial values for mixed states.
You simply have to use a density operator in the function call.

```@example meanfield
u0 = initial_values(me, dm(psi_0))
nothing # hide
```

#### Mapping levels for `NLevelSpace`

The conversion to a numeric representation between [`FockSpace`](@ref) and `FockBasis` is always uniquely defined.
However, there is some freedom of choice when it comes to [`NLevelSpace`](@ref) and the equivalent of `NLevelBasis`, specifically when using symbolic levels.
While it is clear that a symbolic [`Transition`](@ref) operator should map to a numeric `transition`, the choice of which level represents maps to which basis state in the `NLevelBasis` is not fixed.

When using numeric level representations, the [`initial_values`](@ref) and [`to_numeric`](@ref) methods default to using the same numbered basis state:

```@example levelmap
using QuantumCumulants, QuantumOpticsBase
h = NLevelSpace(:TwoLevelAtom, (1, 2))
b = NLevelBasis(2)
s = Transition(h, :s, 1, 2)
@assert to_numeric(s, b) == transition(b, 1, 2)
nothing # hide
```

The order here can be overridden using the `level_map` keyword.
When using symbolic levels, the `level_map` keyword is required.

```@example levelmap2
using QuantumCumulants, QuantumOpticsBase
h = NLevelSpace(:TwoLevelAtom, (:g, :e))
b = NLevelBasis(2)
s = Transition(h, :s, :g, :e)
level_map = Dict(:g => 1, :e => 2)
@assert to_numeric(s, b; level_map=level_map) == transition(b, 1, 2)
nothing # hide
```

#### Numeric averages and conversion

While the examples so far were relatively simple and would have been easy to calculate by hand, things quickly become more difficult whenever product spaces and higher-order products are involved.

Behind the scenes, [`initial_values`](@ref) just uses the [`numeric_average`](@ref) method in order to compute the numeric expectation value for the given operators and states.
This method in turn calls into the numeric conversion [`to_numeric`](@ref) and then uses `QuantumOpticsBase.expect` on the result in order to calculate the respective expectation values for the given state and operators numerically.
Should you need to compute numerical averages from a symbolic one for a given numerical state you can also call [`numeric_average`](@ref) directly.

```@example tonumeric
using QuantumCumulants, QuantumOpticsBase
hfock = FockSpace(:cavity)
hnlevel = NLevelSpace(:ThreeLevelAtom, (:a, :b, :c))
h = hfock ⊗ hnlevel
a = Destroy(h, :a)
s = Transition(h, :s, :a, :c)
levelmap = Dict(
    :a => 3,
    :b => 2,
    :c => 1,
)

bfock = FockBasis(10)
bnlevel = NLevelBasis(3)
psi = coherentstate(bfock, 0.3) ⊗ (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) / sqrt(2)

avg = average(a' * s)
avg_num = numeric_average(avg, psi; level_map=levelmap)
nothing # hide
```

Similarly, you can also just obtain the numerical representation of an operator by directly calling [`to_numeric`](@ref) and a given basis.

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
psi_lazy = LazyKet(b, (coherentstate(bfock, 0.3), (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) / sqrt(2)),)
avg_num_lazy = numeric_average(avg, psi_lazy; level_map=levelmap)
@assert isapprox(avg_num, avg_num_lazy)
```


## [The *q*-number interface](@id interface)

While there are currently only two different Hilbert spaces and two different types of fundamental operators implemented, their implementations are somewhat generic. This means that one can implement custom operator types along with some commutation relations for rewriting. The requirements for that are:

* Custom operator types need to be subtypes of [`QSym`](@ref).
* `Base.:*(::Operator1, ::Operator2)`: A multiplication method that rewrites according to a commutation relation has to be implemented.
* `QuantumCumulants.ismergeable(::Operator1, ::Operator2) = true` is required so pairs of `Operator1` and `Operator2` are detected in longer expressions and rewritten according to their commutation relation.
* Optional: custom Hilbert space type matching the new operators.


#### Example: Harmonic oscillator quadratures

To illustrate, say we would like to implement the quantum harmonic oscillator in terms of the position operator ``x`` and the momentum operator ``p`` rather than the ladder operators. They fulfill the commutation relation

```math
[x,p] = i
```

and we will use it to rewrite occurrences of ``xp \Rightarrow px + i``. For simplicity, we will define them on a [`FockSpace`](@ref) instead of defining a custom Hilbert space as well.

```@example custom-operators
using Latexify # hide
set_default(double_linebreak=true) # hide
using QuantumCumulants

struct Position <: QSym
    hilbert
    name
    aon
    metadata
end
Position(hilbert, name, aon; metadata=QuantumCumulants.source_metadata(:Position, name)) =
    Position(hilbert, name, aon, metadata)

struct Momentum <: QSym
    hilbert
    name
    aon
    metadata
end
Momentum(hilbert, name, aon; metadata=QuantumCumulants.source_metadata(:Momentum, name)) =
    Momentum(hilbert, name, aon, metadata)
nothing # hide
```

Note that any subtype to [`QSym`](@ref) needs to have the four fields shown above, and the 
associated outer constructor. The outer constructor is needed for the interface to 
Symbolics.jl. More fields could be added, but the four shown here are always required. 
Now, for methods we simply need:

```@example custom-operators
using SecondQuantizedAlgebra
SecondQuantizedAlgebra.ismergeable(::Position,::Momentum) = true
Base.:*(x::Position, p::Momentum) = im + p*x
Base.isequal(a::Position, b::Position) = isequal(a.hilbert, b.hilbert) && isequal(a.name, b.name) && isequal(a.aon, b.aon)
Base.isequal(a::Momentum, b::Momentum) = isequal(a.hilbert, b.hilbert) && isequal(a.name, b.name) && isequal(a.aon, b.aon)
```

The `Base.isequal` methods do not compare metadata fields. Note that if your subtypes of 
[`QSym`](@ref) have type parameters, you must also implement a method of `Base.hash` such 
that `isequal(x,y)` implies `hash(x) == hash(y)`.

We can now use our new operator types in expressions and derive equations of motion for them.

```@example custom-operators
h = FockSpace(:oscillator)
x = Position(h,:x,1)
p = Momentum(h,:p,1)

@cnumbers ω m
H = p^2/(2m) + 0.5m*ω^2*x^2

eqs = meanfield([x,p],H)
```
