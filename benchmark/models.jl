# Physical models for the benchmark suite, configs taken from test/systems/ and
# examples/. The pipeline/indexed/corr flags select which feature groups join.

const SQA = QuantumCumulants.SecondQuantizedAlgebra

# Phase-invariance filter, mirrors examples/superradiant_laser_indexed.jl.
φ(x) = 0
φ(::Destroy) = -1
φ(::Create) = 1
φ(x::Transition) = x.i - x.j
function φ(q::SQA.QAdd)
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

function model_jc()
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    @variables Δ g κ γ ν
    H = Δ * a' * a + g * (a' * σ + σ' * a)
    J = [a, σ, σ']
    rates = [κ, γ, ν]
    ops = [a' * a, σ' * σ, a * σ']
    return (;
        name = "jc", ops, H, J, rates, orders = [2, 3], filt = nothing,
        pipeline = true, indexed = false, corr = false,
    )
end

function model_dicke()
    hf = FockSpace(:cavity)
    hs1 = PauliSpace(:spin1)
    hs2 = PauliSpace(:spin2)
    h = hf ⊗ hs1 ⊗ hs2
    a = Destroy(h, :a)
    σ(s, axis) = Pauli(h, Symbol(:σ, s), axis, s + 1)
    @variables Δ_ g κ η
    H = Δ_ * a' * a + g * (a' + a) * (σ(1, 1) + σ(2, 1)) + η * (a' + a)
    J = [a]
    rates = [κ]
    ops = [σ(1, 3), σ(2, 3), σ(1, 3) * σ(2, 3)]
    return (;
        name = "dicke", ops, H, J, rates, orders = [2, 3], filt = nothing,
        pipeline = true, indexed = false, corr = false,
    )
end

function model_multilevel()
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 3)
    h = hf ⊗ ha
    @variables g Δ_2 Δ_3 Ω_2 Ω_3 Δ_c Γ_2 Γ_3 κ
    a = Destroy(h, :a)
    σ(i, j) = Transition(h, :σ, i, j)
    H_atom = -Δ_2 * σ(2, 2) - Δ_3 * σ(3, 3) +
        Ω_2 * (σ(2, 1) + σ(1, 2)) + Ω_3 * (σ(3, 1) + σ(1, 3))
    H_cav = -Δ_c * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
    H = H_atom + H_cav
    J = [a, σ(1, 2), σ(1, 3)]
    rates = [κ, Γ_2, Γ_3]
    ops = [a' * a]
    return (;
        name = "multilevel", ops, H, J, rates, orders = [2, 3], filt = nothing,
        pipeline = true, indexed = false, corr = false,
    )
end

function model_manyatom()
    N = 2  # number of three-level atoms; order-3 over N·NLevel(3) is intractable
    @variables κ g Γ23 Γ13 Γ12 Ω Δc Δ3
    hf = FockSpace(:cavity)
    ha = ⊗([NLevelSpace(Symbol(:atom, i), 3) for i in 1:N]...)
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ(i, j, k) = Transition(h, Symbol("σ_{$k}"), i, j, k + 1)
    H = -Δc * a'a +
        sum(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)) for i in 1:N) +
        sum(Ω * (σ(3, 1, i) + σ(1, 3, i)) for i in 1:N) -
        sum(Δ3 * σ(3, 3, i) for i in 1:N)
    J = [a; [σ(1, 2, i) for i in 1:N]; [σ(1, 3, i) for i in 1:N]; [σ(2, 3, i) for i in 1:N]]
    rates = [κ; [Γ12 for i in 1:N]; [Γ13 for i in 1:N]; [Γ23 for i in 1:N]]
    ops = [a'a, σ(2, 2, 1), σ(3, 3, 1)]
    return (;
        name = "manyatom", ops, H, J, rates, orders = [2], filt = nothing,
        pipeline = true, indexed = false, corr = false,
    )
end

function model_superradiant()
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N Δ κ Γ R ν
    g(i) = IndexedVariable(:g, i)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    H = -Δ * a'a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, j)]
    return (;
        name = "superradiant", ops, H, J, rates, orders = [2, 3],
        filt = phase_invariant, pipeline = true, indexed = true,
        corr = true, op1 = a', op2 = a, steady = true,
    )
end

function model_cavity()
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc κ
    H = ωc * a' * a
    J = [a]
    rates = [κ]
    ops = [a' * a]
    return (;
        name = "cavity", ops, H, J, rates, orders = [2], filt = nothing,
        pipeline = false, indexed = false, corr = true,
        op1 = a, op2 = a', steady = false,
    )
end

models() = [
    model_jc(),
    model_dicke(),
    model_multilevel(),
    model_manyatom(),
    model_superradiant(),
    model_cavity(),
]
