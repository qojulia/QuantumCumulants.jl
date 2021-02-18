# Implementation

Let's take a closer look at each step involved from defining a system to arriving at a numerical solution of the underlying time dynamics.


## Hilbert spaces

The first step in treating a system with **Qumulants.jl** is to specify the Hilbert space on which the system is defined. There are two types of Hilbert spaces implemented, namely [`FockSpace`](@ref) and [`NLevelSpace`](@ref). The first describes systems whose operators follow the fundamental bosonic commutation relations (such as the quantum harmonic oscillator), whereas the latter describes systems consisting of a number of discrete energy levels (such as atoms).

A [`FockSpace`](@ref) simply needs a name in order to be defined:

```@example hilbert-space
using Qumulants # hide
hf = FockSpace(:fock1)
nothing # hide
```

[`NLevelSpace`](@ref) requires a name as well as labels for the energy levels. For example

```@example hilbert-space
h_atom = NLevelSpace(:atom, (:g,:e))
nothing # hide
```

defines an [`NLevelSpace`](@ref) with the name `:atom` and the two levels labeled by `:g` and `:e`, respectively. Note that the levels can be labeled by (almost) anything. For example, `NLevelSpace(:two_level, (1,2))` would define a Hilbert space describing a system with the two discrete energy lavels labeled by `1` and `2`. Specifically for numbers, there is also the short-hand method to write `NLevelSpace(:five_level, 5)` which creates a system with levels `1:5`. Note that by default the first level in the list of all levels is designated as the ground state. This can be changed by specifying the ground state explicitly as a third argument to [`NLevelSpace`](@ref). The ground state projector will be eliminated during simplification (see below).

Composite systems are generally described by a [`ProductSpace`](@ref), i.e. a Hilbert space that consists of multiple subspaces. Each subspace is either a [`FockSpace`](@ref) or an [`NLevelSpace`](@ref). They can be created using the [`tensor`](@ref) function or the unicode symbol [`⊗`](@ref) [\otimes]. For example

```@example hilbert-space
h_prod1 = tensor(hf, h_atom)
h_prod2 = tensor(h_prod1, NLevelSpace(:three_level, 3))
nothing # hide
```

creates two product spaces. The first, `h_prod1`, consists of the previously defined `FockSpace(:fock1)` and `NLevelSpace(:atom, (:g,:e))`. The second one, `h_prod2`, adds in another `NLevelSpace(:three_level, 3)`. In principle arbitrarily many systems can be combined this way.


## Operators (a.k.a. *q*-numbers)

Once the Hilbert space of the system has been defined, we can proceed by defining operators, or *q*-numbers, on them. They are the fundamental building blocks of symbolic expressions in **Qumulants.jl**. Again, there are essentially two kinds of operators implemented: the quantum harmonic destruction operator [`Destroy`](@ref) which acts on a [`FockSpace`](@ref), as well as a [`Transition`](@ref) operator which describes a transition between any two levels on an [`NLevelSpace`](@ref). Of course, these operators can only be defined on the corresponding Hilbert spaces.

Here are a few examples:

```@example operators
using Qumulants # hide
hf = FockSpace(:fock)
a = Destroy(hf, :a)

h_atom = NLevelSpace(:atom,(:g,:e))
σge = Transition(h_atom, :σ, :g, :e)
σ = Transition(h_atom, :σ)
@assert isequal(σge, σ(:g,:e)) # true
nothing # hide
```

As you can see, the destruction operator [`Destroy`](@ref) is created on a [`FockSpace`](@ref) and given a name. The transition operator, however, additionally requires you to specify the levels between which it describes the transition. Defining a transition without levels specified creates a callable instance which needs to be called with valid level labels before one can actually use it in any algebraic expressions. Note that in Bra-Ket notation, the transition operator `Transition(h, i, j)` is simply `|i><j|`.

These fundamental operators are all of the type [`QSym`](@ref), which are the basic symbolic building blocks for the noncommutative algebra used in **Qumulants.jl**. They can be combined with standard algebraic functions in expression trees, which are implemented as [`QTerm`](@ref).

```@example operators
ex_fock = 0.1*a'*a
ex_trans = im*(σ(:g,:e) - σ(:e,:g))
nothing # hide
```

Note that only operators that are defined on the same Hilbert space can be algebraically combined.

In composite systems, we also need to specify on which subsystem the respective operator acts. This information is important as operators acting on different subsystems commute with one another, but operators acting on the same one do not. When multiplying together operators in a composite systems, they are automatically ordered according to the order of Hilbert spaces.

```@example operators
h_prod = FockSpace(:fock1) ⊗ FockSpace(:fock2)
a = Destroy(h_prod,:a,1)
b = Destroy(h_prod,:b,2)
a*b # a*b
b*a # b*a
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

Commutative numbers (*c*-numbers) are represented by `SymbolicUtils.Sym` from the [**SymbolicUtils.jl**](https://github.com/JuliaSymbolics/SymbolicUtils.jl) and a custom subtype to `Number` called [`CNumber`](@ref). They are generally assumed to be complex numbers and can be defined with the [`cnumbers`](@ref) function or the corresponding macro [`@cnumbers`](@ref). You can use them together with *q*-numbers to build symbolic expressions describing the Hamiltonian, e.g.

```@example c-numbers
using Qumulants # hide
h = FockSpace(:fock)
@cnumbers ω η
@qnumbers a::Destroy(h)
H = ω*a'*a + η*(a + a')
nothing # hide
```


## Deriving Heisenberg equations and simplification

The equations of motion of *q*-numbers are determined by evaluating commutators. This can be done by using fundamental commutation relations in the term rewriting rules.

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
\sigma^{ij}\sigma^{kl} = \delta_{jk}\sigma^{il},
```

which is implemented as rewriting rule just so. Additionally, we use the fact that in a system with levels ``\{1,...,n\}``

```math
\sum_{j=1}^n \sigma^{jj} = 1
```

in order to eliminate the projector on the ground state. This reduces the amount of equations required for each [`NLevelSpace`](@ref) by 1. Note that, as mentioned before, the ground state is by default chosen to be the first (but this can be cahnged). Hence, the default rewriting rule to eliminate the ground-state projector is

```math
\sigma^{11} ~\Rightarrow~ 1 - \sum_{j=2}^n \sigma^{jj}.
```

These rules are applied automatically when deriving Heisenberg equations for *q*-numbers. If you manually want to apply them to an expression you can use [`qsimplify`](@ref):

```@example heisenberg
using Qumulants # hide
h = FockSpace(:fock)
@qnumbers a::Destroy(h)
qsimplify(a*a') # return a'*a + 1
nothing # hide
```

In order to derive equations of motion, you need to specify a Hamiltonian and the operator (or a list of operators) of which you want to derive the Heisenberg equations and pass them to [`heisenberg`](@ref).

```@example heisenberg
using Latexify # hide
set_default(double_linebreak=true) # hide
@cnumbers ω η
H = ω*a'*a + η*(a + a') # Driven cavity Hamiltonian
he = heisenberg([a, a'*a], H)
```

To add decay to the system, you can pass an additional list of operators corresponding to the collapse operators describing the respective decay. For example, `heisenberg(a, H, [a]; rates=[κ])` would derive the equations of a cavity that is also subject to decay at a rates `κ`.

The equations resulting from the call to [`heisenberg`](@ref) are stored as an instance of [`HeisenbergEquation`](@ref), which stores the left-hand-side and the right-hand-side of the equations together with additional information such as the system Hamiltonian.


## Cumulant expansion

Averaging (using [`average`](@ref)) and the [`cumulant_expansion`](@ref) are essential to convert the system of *q*-number equations to *c*-number equations. Averaging alone converts any operator product to a *c*-number, yet you will not arrive at a closed set of equations without truncating at a specific order. An average is stored as a symbolic expression. Specifically, the average of an operator `op` is internally represented by `SymbolicUtils.Term{AvgSym}(average, [op])`.

The order of an average is given by the number of constituents in the product. For example

```@example cumulant
using Qumulants # hide
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
average(a'*a, 1) # short-hand for cumulants_expansion(average(a'*a), 1)
nothing #hide
```

When applying this to a system of equations, you obtain a set of *c*-number differential equations.

```@example heisenberg
he_avg = average(he,2)
```

Before you can actually solve the system of equations, you need to ensure that it is complete, i.e. there are no averages missing. This can be checked with [`find_missing`](@ref). Alternatively, you can automatically complete a system of equations using the [`complete`](@ref) function which will internally use [`find_missing`](@ref) to look for missing averages and derive equations for those.


## Numerical solution

Finally, in order to actually solve a system of equations, we need to generate a function that can be used in the [**OrdinaryDiffEq.jl**](https://github.com/SciML/OrdinaryDiffEq.jl) package. This is done using the [`build_ode`](@ref) function, which generates an `Expr` that can be evaluated, or the [`generate_ode`](@ref) which calls `Meta.eval` on the result of [`build_ode`](@ref).

```@example heisenberg
ps = (ω, η)
meta_f = build_ode(he_avg, ps)
using MacroTools; MacroTools.striplines(meta_f)
```

Note that in addition to the system, we also need to pass a list of symbolic parameters. The order in that list is important as it will correspond to the order in which parameters need to be passed. Similarly, the solution vector will have entries in the order of the left-hand-side of the equations. The resulting function is of the form `f(du,u,p,t)`, which is precisely the format required in [**OrdinaryDiffEq.jl**](https://github.com/SciML/OrdinaryDiffEq.jl). We can solve it as follows

```@example heisenberg
using OrdinaryDiffEq
f = generate_ode(he_avg, ps) # Meta.eval(meta_f)
p0 = (1.0, 0.1) # (ω, η) in that order
u0 = zeros(ComplexF64, length(he_avg))
prob = ODEProblem(f,u0,(0.0,1.0),p0)
sol = solve(prob, RK4())
nothing # hide
```

Now, the state of the system at each time-step is stored in `sol.u`. As mentioned above, the order of the solution is the same as in the left-hand-side of `he_avg`. For example, `sol.u[end][1]` corresponds to ``\langle a \rangle`` at the final time. In large systems it might be tedious to keep track of the order of the expectation values, especially when using [`complete`](@ref). For convenience, you can use [`get_solution`](@ref) in such a case.
