```@meta
CollapsedDocStrings = true
```

# [API](@id API)

```@contents
Pages = ["API.md"]
Depth = 2:3
```

## [Hilbert Spaces](@id API: Hilbert Spaces)

```@docs
SecondQuantizedAlgebra.HilbertSpace
```

```@docs
ProductSpace
```

```@docs
FockSpace
```

```@docs
NLevelSpace
```

```@docs
SecondQuantizedAlgebra.CollectiveNLevelSpace
```

```@docs
PauliSpace
```

```@docs
SpinSpace
```

```@docs
PhaseSpace
```

```@docs
⊗
```

```@docs
tensor
```

## [q-Numbers](@id API: Operators)

```@docs
SecondQuantizedAlgebra.QField
```

```@docs
SecondQuantizedAlgebra.QSym
```

```@docs
SecondQuantizedAlgebra.QAdd
```

```@docs
SecondQuantizedAlgebra.QTerm
```

```@docs
SecondQuantizedAlgebra.QTermDict
```

```@docs
SecondQuantizedAlgebra.Op
```

```@docs
SecondQuantizedAlgebra.OpKind
```

```@docs
SecondQuantizedAlgebra.get_prefactor
SecondQuantizedAlgebra.get_operators
SecondQuantizedAlgebra.get_variables
SecondQuantizedAlgebra.sorted_arguments
SecondQuantizedAlgebra.constraint_pairs
```

```@docs
SecondQuantizedAlgebra.normal_order
SecondQuantizedAlgebra.qadjoint
SecondQuantizedAlgebra.inner_adjoint
SecondQuantizedAlgebra.symmetric_to_normal
SecondQuantizedAlgebra.normal_to_symmetric
SecondQuantizedAlgebra.expand_completeness
SecondQuantizedAlgebra.anticommutator
SecondQuantizedAlgebra.is_average
SecondQuantizedAlgebra.is_indexed_sum
SecondQuantizedAlgebra.undo_average
SecondQuantizedAlgebra.dagger(::SecondQuantizedAlgebra.QField)
SymbolicUtils.substitute
SymbolicUtils.simplify
SymbolicUtils.expand
```

```@docs
@qnumbers
```

```@docs
Destroy
```

```@docs
Create
```

```@docs
Transition
```

```@docs
SecondQuantizedAlgebra.CollectiveTransition
```

```@docs
Pauli
```

```@docs
Spin
```

```@docs
Position
```

```@docs
Momentum
```

```@docs
SecondQuantizedAlgebra.optype
SecondQuantizedAlgebra.operator_name
SecondQuantizedAlgebra.operator_index
SecondQuantizedAlgebra.is_destroy
SecondQuantizedAlgebra.is_create
SecondQuantizedAlgebra.is_transition
SecondQuantizedAlgebra.is_collective_transition
SecondQuantizedAlgebra.is_pauli
SecondQuantizedAlgebra.is_spin
SecondQuantizedAlgebra.is_position
SecondQuantizedAlgebra.is_momentum
```

## [Unitary Transformations](@id API: Unitary)

```@docs
SecondQuantizedAlgebra.UnitaryTransform
SecondQuantizedAlgebra.conjugate
SecondQuantizedAlgebra.transform
SecondQuantizedAlgebra.Displace
SecondQuantizedAlgebra.DisplacementFrame
SecondQuantizedAlgebra.Rotation
SecondQuantizedAlgebra.RotatingFrame
SecondQuantizedAlgebra.Squeeze
SecondQuantizedAlgebra.Bogoliubov
SecondQuantizedAlgebra.gauge_term
SecondQuantizedAlgebra.generators
```

## [Mean field](@id API: Meanfield)

```@docs
meanfield
```

```@docs
commutator
```

```@docs
acts_on
```

```@docs
AbstractMeanfieldEquations
```

```@docs
MeanfieldEquations
```

## [Average](@id API: Average)

```@docs
average
```

```@docs
cumulant_expansion
```

```@docs
cumulant
```

```@docs
get_order
```

## [Introspection](@id API: Introspection)

```@docs
states
```

```@docs
operators
```

```@docs
moments
```

```@docs
moment_variable_map
```

```@docs
closure_report
```

```@docs
noise_channels
```

## [Correlation functions](@id API: correlation)

```@docs
CorrelationFunction
```

```@docs
ModelingToolkitBase.System
```

```@docs
Spectrum
```

```@docs
correlation_u0
```

```@docs
correlation_p0
```

## [Symbolic Summations](@id API: Sums)

```@docs
Index
```

```@docs
IndexedOperator
```

```@docs
IndexedVariable
```

```@docs
DoubleIndexedVariable
```

```@docs
Σ
```

```@docs
change_index
```

```@docs
get_indices
SecondQuantizedAlgebra.has_index
SecondQuantizedAlgebra.assume_distinct_index
SecondQuantizedAlgebra.index_slot
SecondQuantizedAlgebra.index_range
SecondQuantizedAlgebra.index_name
SecondQuantizedAlgebra.index_sym
SecondQuantizedAlgebra.get_sum_indices
SecondQuantizedAlgebra.get_sum_non_equal
SecondQuantizedAlgebra.has_sum_metadata
```

```@docs
evaluate
```

```@docs
scale
```

```@docs
scale!
```

## [Measurement Backaction](@id API: Measurement Backaction)

```@docs
NoiseMeanfieldEquations
```

```@docs
QuantumCumulants.translate_W_to_Y
```

```@docs
EvolutionDirection
```

```@docs
Forward
```

```@docs
Backward
```

## [Utility functions](@id API: Utils)

```@docs
find_missing
```

```@docs
find_operators
```

```@docs
complete
```

```@docs
complete!
```

```@docs
unique_up_to_adjoint
SecondQuantizedAlgebra.unique_up_to_adjoint!
```

```@docs
QuantumCumulants.simplify!
```

```@docs
fundamental_operators
```

```@docs
to_numeric
```

```@docs
numeric_average
```

```@docs
initial_values
```

```@docs
get_solution
```

```@docs
parameter_map
```

```@docs
modify_equations
```

```@docs
modify_equations!
```

```@docs
substitute!
```

```@docs
SecondQuantizedAlgebra.make_time_dependent
```

## [Direct ODE backend](@id API: Direct ODE backend)

For completed deterministic moment hierarchies, the direct backend bypasses construction and
compilation of a ModelingToolkit `System`. The equations are lowered once to a structured
polynomial representation and executed by a compact numerical kernel with a concrete
parameter payload.

A typical solve is:

```julia
using QuantumCumulants
using SciMLBase: ODEProblem
using OrdinaryDiffEqTsit5: Tsit5, solve

closed = complete(eqs)
u0 = zeros(ComplexF64, length(closed.states))
ps = Dict(g => 1.0, κ => 0.2)

prob = ODEProblem(
    closed,
    u0,
    (0.0, 10.0),
    ps;
    backend = KernelBackend(),
)
sol = solve(prob, Tsit5())
```

The backend is explicitly opt-in. It accepts polynomial deterministic drift with
state-independent coefficients. Time-dependent coefficients or unsupported state-dependent
functions should continue through the general `System(eqs)` path.

Parameter values can be changed without relowering the equations or rebuilding the numerical
kernel:

```julia
update_parameters!(prob, Dict(g => 1.2))
sol2 = solve(prob, Tsit5())
```

Updates may be partial: unspecified parameters retain their current values. For independent
parameter sweeps, copy or remake the problem rather than mutating one shared problem from
multiple tasks.

```@docs
KernelBackend
update_parameters!
```

## [Numeric backends](@id API: Numeric backends)

Symbolic-to-numeric conversion ([`to_numeric`](@ref), [`numeric_average`](@ref),
[`initial_values`](@ref)) runs through SecondQuantizedAlgebra's pluggable backend
interface. Loading a backend (`using QuantumOpticsBase` or `using QuantumToolbox`) selects
the concrete implementation; the conversion details are documented in SQA's
[numeric conversion section](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/implementation/#Numeric-conversion).
The re-exported types and hooks below are the interface for selecting and implementing a
backend.

```@docs
NumericBackend
QuantumOpticsBackend
QuantumToolboxBackend
numeric_backend
numeric_basis
numeric_subbasis
numeric_operator
numeric_identity
numeric_embed
numeric_num_subsystems
numeric_assemble
numeric_assemble_td
numeric_materialize
numeric_expect
```
