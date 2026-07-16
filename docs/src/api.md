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
SecondQuantizedAlgebra.QSym
```

```@docs
SecondQuantizedAlgebra.QAdd
```

```@docs
SecondQuantizedAlgebra.QTerm
```

```@docs
SecondQuantizedAlgebra.Op
```

```@docs
SecondQuantizedAlgebra.prefactor
SecondQuantizedAlgebra.operators
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
SecondQuantizedAlgebra.is_destroy
SecondQuantizedAlgebra.is_create
SecondQuantizedAlgebra.is_transition
SecondQuantizedAlgebra.is_pauli
SecondQuantizedAlgebra.is_spin
SecondQuantizedAlgebra.is_position
SecondQuantizedAlgebra.is_momentum
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

## [RHS backends](@id API: RHS backends)

See [Solving the equations directly (RHS backends)](@ref) for the guide.

```@docs
RHSBackend
```

```@docs
AutoBackend
```

```@docs
KernelBackend
```

```@docs
ShardedBackend
```

```@docs
update_parameters!
```

```@docs
SciMLBase.ODEFunction(::MeanfieldEquations, ::Any)
```

```@docs
SciMLBase.ODEProblem(::MeanfieldEquations, ::Any, ::Any, ::Any)
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
unique_ops
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
SecondQuantizedAlgebra.make_time_dependent
```
