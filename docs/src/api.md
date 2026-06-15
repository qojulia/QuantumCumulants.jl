# [API](@id API)

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
