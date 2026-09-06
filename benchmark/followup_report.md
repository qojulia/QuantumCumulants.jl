# Structured backend follow-up

These measurements were run against QuantumCumulants checkpoint `2216582` with
Julia 1.12.7, 12 Julia threads, `BenchmarkTools`, `OrdinaryDiffEqTsit5/Tsit5`,
and `RuntimeGeneratedFunctions` on the N=6 transverse-field Ising fixture. Each cold result used a
fresh Julia process or a fresh persistent Julia worker. Package loading and
precompilation were outside the measured pipeline. Compilation triggered by a
stage is included in that stage. Warm RHS numbers are separate steady-state
measurements taken only after the cold construction and first call completed.

The measured cold pipeline is:

```text
meanfield -> complete -> [MomentIR lowering or MTK System] -> evaluator/kernel
-> parameter payload -> ODEFunction/ODEProblem -> first Tsit5 solution
```

The direct and generated paths use `MomentIR` with shared monomial IDs. The
generated experiment packs consecutive units by operation count: monomial cost
is the number of factors and row cost is the number of sparse terms. Generated
units are held in `Vector{Function}` and called through a `@noinline` wrapper,
so the dispatcher cannot fold all units into one giant method. The native
experiment is confined to this benchmark directory.

## Cold path results

Times are seconds. `after complete` includes every stage from the next row of
the pipeline through the first solution. It excludes `meanfield` and
`complete`, which are reported separately because they are common equation
construction stages.

| equations | path | lowering/system | evaluator or compile | parameter payload | first solution | total after complete |
|---:|---|---:|---:|---:|---:|---:|
| 693 | compact | MomentIR 2.648 | kernel + ODEFunction 0.054 | 2.135 | 2.690 | 7.526 |
| 693 | MTK | System 3.301 + compile 2.957 | ODEProblem 13.235 | 0.050 | 49.774 | 69.317 |
| 693 | generated, target 64 | setup 5.064 | first call compile 34.002 | included | 1.957 | 41.023 |
| 693 | generated, target 128 | setup 5.116 | first call compile 31.364 | included | 1.965 | 38.446 |
| 693 | generated, target 256 | setup 5.103 | first call compile 33.380 | included | 1.967 | 40.451 |
| 693 | generated rows only, target 256 | setup 5.242 | first call compile 15.143 | included | 1.975 | 22.360 |
| 1908 | compact | MomentIR 12.729 | kernel + ODEFunction 0.057 | 2.099 | 2.551 | 17.435 |
| 1908 | generated rows only, target 256 | setup 15.611 | first call compile 100.537 | included | 1.954 | 118.101 |

The order-3 common stages were `meanfield=13.106 s` and `complete=1.197 s`
for compact, and `meanfield=12.353 s` and `complete=1.259 s` for MTK. The
order-4 compact common stages were `meanfield=13.245 s` and
`complete=3.694 s`. The equation counts were 693 and 1,908; the corresponding
monomial counts were 12,341 and 66,451.

## Compact RHS profile

These are median steady-state timings from 100 one-call samples after warmup.
Allocation counts and bytes are from the minimum-allocation sample.

| equations | scratch lookup | monomial pass | sparse accumulation | complete RHS |
|---:|---:|---:|---:|---:|
| 693 | 0.110 μs, 1 alloc / 16 B | 9.82 μs, 0 / 0 B | 12.06 μs, 0 / 0 B | 23.3 μs, 2 / 64 B |
| 1908 | 0.110 μs, 1 alloc / 16 B | 60.9 μs, 0 / 0 B | 85.5 μs, 0 / 0 B | 142.7 μs, 2 / 64 B |

Sparse accumulation is the larger of the two numeric passes at both sizes, so
the focused native experiment generated rows only as well as the full
monomial-plus-row variant. The row-only generated RHS was 74.6 μs with 83
allocations / 21,336 B at order 3 and 776 μs with 650 allocations / 113,304 B
at order 4. It did not improve the compact evaluator.

The one-second integrations from nonzero complex initial states took 0.859 s
for 693 equations and 0.870 s for 1,908 equations, each saving 101 points.
These are solve timings, not RHS-only timings, and use the same Tsit5 setup as
the short first-solution measurements.

## Numerical checks

For every generated order-3 target and the order-4 row-only run:

- RHS values matched compact at three deterministic nonzero complex states,
  with reported maximum error `0.0`.
- After changing one parameter, the maximum RHS error remained `0.0`.
- An 11-point Tsit5 trajectory from a nonzero complex state had maximum error
  `0.0` against compact.

The equality is expected here because both evaluators use the same numeric
parameter payload and operation order. It checks the generated path against
compact at several states and after a parameter update; it does not establish
agreement with an independent symbolic or physical reference model.

## MTK order 4 limit

The order-4 MTK run was given a separate Julia worker thread and durable output.
It reached 1,908 completed equations, then remained in MTK system construction
for 600 seconds. RSS stayed below the explicit 6 GiB ceiling, peaking around
5.24 GiB. The worker was terminated at the bound. Therefore there is no
order-4 MTK first-solution number to report. The durable record is
`results/followup-mtk-order4-thread.log`; the earlier 300-second MCP attempt is
also retained in `results/followup-mtk-order4.log`.

## Findings

The compact backend is the useful production path in the measured range. Its
first solution was 2.69 s at 693 equations and 2.55 s at 1,908 equations, with
warm RHS costs of 23.3 μs and 142.7 μs and no allocations in either isolated
numeric pass. At order 3, MTK took 69.3 s from completed equations to its first
solution and its warm RHS was 32.1 μs.

The corrected generated design preserved numerical behavior but added large
compile and allocation costs. At order 3, operation targets 64, 128, and 256
gave 38.4--41.0 s to the first solution and remained slower than compact in
warm RHS. At order 4, even the focused row-only design took 118.1 s cold,
100.5 s in first-call compilation, and 0.776 ms per warm RHS. The measured
native design is therefore discarded for the main backend. This verdict applies
to the operation-bounded, noinline-dispatch design tested here.

The order-4 MomentIR lowering allocated 44.38 GB cumulatively in 12.73 s.
That is a profiling lead rather than a production optimization target for this
checkpoint. Peak RSS was not inferred from cumulative allocation counters.

## Reproduction

Run each command in a fresh Julia process from the repository root:

```text
QC_FOLLOWUP_MODE=compact QC_FOLLOWUP_ORDER=3 \
  julia --project=benchmark benchmark/followup_backend.jl
QC_FOLLOWUP_MODE=compact QC_FOLLOWUP_ORDER=4 \
  julia --project=benchmark benchmark/followup_backend.jl
QC_FOLLOWUP_MODE=mtk QC_FOLLOWUP_ORDER=3 \
  julia --project=benchmark benchmark/followup_backend.jl
QC_FOLLOWUP_MODE=generated QC_FOLLOWUP_ORDER=3 QC_GEN_TARGET=256 \
  julia --project=benchmark benchmark/followup_backend.jl
QC_FOLLOWUP_MODE=generated QC_FOLLOWUP_ORDER=4 QC_GEN_STAGE=rows QC_GEN_TARGET=256 \
  julia --project=benchmark benchmark/followup_backend.jl
```

The runner writes progress to `benchmark/results/` through
`QC_FOLLOWUP_RESULTS`. The generated order-3 target sweep is summarized in
`results/generated-order3-targets.tsv`.
