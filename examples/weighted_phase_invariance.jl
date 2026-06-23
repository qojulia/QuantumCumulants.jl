# # Weighted Phase Invariance

# Many quantum optical systems possess a $U(1)$ phase symmetry: the Hamiltonian is invariant under $a \rightarrow a\, e^{i\theta}$ for every mode, which implies that every average carrying a net phase vanishes. The [Superradiant Laser](@ref) example exploits this to discard all phase-carrying averages with a `filter_func`, dramatically shrinking the set of equations.

# Sometimes the conserved quantity is not the bare excitation number but a *weighted* excitation number. Consider three bosonic modes $a$, $b$ and $f$ coupled by a non-degenerate down-conversion process in which one $a$ photon splits into one $b$ photon and one $f$ photon,
# ```math
# H = \Delta_a a^\dagger a + \Delta_b b^\dagger b + \Delta_f f^\dagger f + g\,(a^\dagger b f + a\, b^\dagger f^\dagger).
# ```
# Energy conservation of the interaction ties the mode frequencies together, and the dynamics conserves the weighted excitation number
# ```math
# 2\, a^\dagger a + b^\dagger b + f^\dagger f ,
# ```
# i.e. mode $a$ carries weight $2$ while $b$ and $f$ carry weight $1$. The symmetry is now $a \rightarrow a\, e^{2 i\theta}$, $b \rightarrow b\, e^{i\theta}$, $f \rightarrow f\, e^{i\theta}$, and an average survives only if its weights sum to zero. This is exactly the question raised in [issue #237](https://github.com/qojulia/QuantumCumulants.jl/issues/237): how to write the filter when the modes carry different weights.

# We start by loading the packages.

using QuantumCumulants
using ModelingToolkitBase, OrdinaryDiffEqTsit5
using Plots

# Each mode lives on its own Fock space. We keep track of which subspace an operator acts on through its space index (the second argument of the operator constructor), since the weight is attached to the mode, not to the operator type.

ha = FockSpace(:a)
hb = FockSpace(:b)
hf = FockSpace(:f)
h = ha ⊗ hb ⊗ hf

@qnumbers a::Destroy(h, 1) b::Destroy(h, 2) f::Destroy(h, 3)
@variables Δa Δb Δf g κa κb κf

H = Δa * a' * a + Δb * b' * b + Δf * f' * f + g * (a' * b * f + a * b' * f')

J = [a, b, f] # photon loss from each mode
rates = [κa, κb, κf]
nothing # hide

# We now build the filter. The phase function $\varphi$ assigns to each operator its weight under the $U(1)$ rotation. The only change compared to the unweighted case is that an annihilation or creation operator contributes $\mp w$, where $w$ is the weight of the mode it acts on. We read the mode off the operator's `space_index` and look the weight up in a dictionary.

import QuantumCumulants.SecondQuantizedAlgebra as SQA

weights = Dict(1 => 2, 2 => 1, 3 => 1) # space_index => U(1) weight

φ(x) = 0
function φ(op::SQA.Op)
    SQA.is_destroy(op) && return -weights[op.space_index]
    SQA.is_create(op) && return weights[op.space_index]
    return 0
end
function φ(q::SQA.QAdd) # walk operator-product expression
    for (term, _) in q.arguments
        p = 0
        for op in term.ops
            p += φ(op)
        end
        return p
    end
    return 0
end
function φ(avg)
    SQA.is_average(avg) || return 0
    return φ(SQA.undo_average(avg))
end
phase_invariant(x) = iszero(φ(x))

# Only the weights' ratios matter, so `Dict(1 => 4, 2 => 2, 3 => 2)` would give the same filter. We can check that the interaction term and the relevant correlators are kept while single-mode coherences are discarded.

φ(average(a' * b * f)), φ(average(a)), φ(average(a' * b))

# The interaction $a^\dagger b f$ involves three operators, so a closed and non-trivial description requires a third-order cumulant expansion. We derive the equations for the mode populations and complete the system with the filter.

ops = [a' * a, b' * b, f' * f]
eqs = meanfield(ops, H, J; rates = rates, order = 3)
eqs_c = complete(eqs; filter_func = phase_invariant)
nothing # hide

# The filter keeps only the weighted-phase-invariant averages. Without it the completion is far larger.

length(eqs_c), length(complete(eqs; order = 3))

# Finally we evolve the system. We start with four photons in mode $a$ and watch them convert into pairs of $b$ and $f$ photons before everything decays.

@named sys = System(eqs_c)
sys = complete(sys)

u0 = zeros(ComplexF64, length(eqs_c))
u0[1] = 4.0 # ⟨a'a⟩(0) = 4, all other averages zero
ps = [Δa => 0.0, Δb => 0.0, Δf => 0.0, g => 1.0, κa => 0.1, κb => 0.2, κf => 0.2]

prob = ODEProblem(sys, merge(initial_values(eqs_c, u0), Dict(ps)), (0.0, 30.0))
sol = solve(prob, Tsit5())
nothing # hide

# We extract the populations and plot them. The conversion of one $a$ photon into a $b$ and an $f$ photon keeps $\langle b^\dagger b \rangle$ and $\langle f^\dagger f \rangle$ equal at all times.

na = real.(get_solution(sol, a' * a, eqs_c).(sol.t))
nb = real.(get_solution(sol, b' * b, eqs_c).(sol.t))
nf = real.(get_solution(sol, f' * f, eqs_c).(sol.t))

plot(sol.t, na; label = "⟨a'a⟩", xlabel = "t", ylabel = "photon number")
plot!(sol.t, nb; label = "⟨b'b⟩")
plot!(sol.t, nf; label = "⟨f'f⟩", linestyle = :dash)
