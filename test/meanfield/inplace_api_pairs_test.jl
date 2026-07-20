using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
using ModelingToolkitBase: System, equations, unknowns
using Test

# Every `foo` / `foo!` pair must be semantically identical: `foo` returns what
# `foo!` produces on a fresh copy, and `foo` must not touch its input. These
# tests pin both properties for `complete`, `scale`, `simplify`, and
# `modify_equations`, on deterministic and noise equations, and confirm the
# downstream `System` agrees.

# ---- semantic comparison helpers --------------------------------------------

_rhs_equal(p, q) = isequal(p.lhs, q.lhs) && _is_zero(p.rhs - q.rhs)

"""
    assert_pair_identical(A, B)

Assert two equation sets are semantically identical: equations (LHS and RHS),
states, operators, per-subspace treatments, the noise drift (when both are
noise sets), and the equations of the assembled `System`.
"""
function assert_pair_identical(A, B)
    @test typeof(A) == typeof(B)
    @test length(A.equations) == length(B.equations)
    @test all(_rhs_equal(p, q) for (p, q) in zip(A.equations, B.equations))
    @test length(A.states) == length(B.states)
    @test all(isequal(p, q) for (p, q) in zip(A.states, B.states))
    @test length(A.operators) == length(B.operators)
    @test all(isequal(p, q) for (p, q) in zip(A.operators, B.operators))
    @test isequal(A.graph.treatments, B.graph.treatments)
    if A isa NoiseMeanfieldEquations && B isa NoiseMeanfieldEquations
        @test length(A.noise_equations) == length(B.noise_equations)
        @test all(_rhs_equal(p, q) for (p, q) in zip(A.noise_equations, B.noise_equations))
    end
    sysA = System(A; name = :s)
    sysB = System(B; name = :s)
    @test length(equations(sysA)) == length(equations(sysB))
    @test length(unknowns(sysA)) == length(unknowns(sysB))
    @test all(isequal(x.lhs, y.lhs) for (x, y) in zip(equations(sysA), equations(sysB)))
    return all(_is_zero(x.rhs - y.rhs) for (x, y) in zip(equations(sysA), equations(sysB)))
end

"""
    assert_unchanged(A, snapshot)

Assert the input `A` was not mutated by a non-mutating call: it still matches a
freshly built reference `snapshot`.
"""
function assert_unchanged(A, snapshot)
    @test length(A.equations) == length(snapshot.equations)
    @test all(_rhs_equal(p, q) for (p, q) in zip(A.equations, snapshot.equations))
    return @test isequal(A.graph.treatments, snapshot.graph.treatments)
end

# ---- fixtures ---------------------------------------------------------------

# anharmonic cavity: incomplete at order 2, so `complete` has work to do
function _det_eqs()
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ U
    H = ω * a' * a + U * a' * a' * a * a
    return meanfield([a], H, [a]; rates = [κ], order = 2)
end

function _noise_eqs()
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ U η
    H = ω * a' * a + U * a' * a' * a * a
    return meanfield([a, a' * a], H, [a]; rates = [κ], efficiencies = [η], order = 2)
end

# indexed Tavis-Cummings: needs `scale` to reduce the symmetric atom subspace
function _indexed_eqs(; efficiencies = nothing)
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @variables Δ g γ N
    i = Index(h, :i, N, ha)
    c = Destroy(h, :c)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)
    H = Δ * c' * c + g * Σ(σ(2, 1, i) * c + σ(1, 2, i) * c', i)
    rates = efficiencies === nothing ? [γ] : [γ]
    return meanfield(c, H, [c]; rates = rates, efficiencies = efficiencies, order = 2)
end

# ---- complete / complete! ---------------------------------------------------

@testset "complete vs complete! (deterministic)" begin
    nonmut = complete(_det_eqs())
    input = _det_eqs()
    mut = complete!(input)
    assert_pair_identical(nonmut, mut)
end

@testset "complete leaves its input untouched" begin
    input = _det_eqs()
    snapshot = _det_eqs()
    complete(input)
    assert_unchanged(input, snapshot)
end

@testset "complete vs complete! (noise)" begin
    nonmut = complete(_noise_eqs())
    @test nonmut isa NoiseMeanfieldEquations
    input = _noise_eqs()
    mut = complete!(input)
    assert_pair_identical(nonmut, mut)
end

# ---- scale / scale! ---------------------------------------------------------

@testset "scale vs scale! (deterministic, indexed)" begin
    nonmut = scale(complete(_indexed_eqs()))
    input = complete(_indexed_eqs())
    mut = scale!(input)
    assert_pair_identical(nonmut, mut)
end

@testset "scale leaves its input untouched" begin
    input = complete(_indexed_eqs())
    snapshot = complete(_indexed_eqs())
    scale(input)
    assert_unchanged(input, snapshot)
end

@testset "scale vs scale! (noise, indexed)" begin
    nonmut = scale(complete(_indexed_eqs(; efficiencies = [1 // 2])))
    @test nonmut isa NoiseMeanfieldEquations
    input = complete(_indexed_eqs(; efficiencies = [1 // 2]))
    mut = scale!(input)
    assert_pair_identical(nonmut, mut)
end

# ---- simplify / simplify! ---------------------------------------------------

@testset "simplify vs simplify! (deterministic)" begin
    base() = complete(_det_eqs())
    nonmut = SymbolicUtils.simplify(base())
    input = base()
    mut = simplify!(input)
    assert_pair_identical(nonmut, mut)
end

@testset "simplify leaves its input untouched" begin
    input = complete(_det_eqs())
    snapshot = complete(_det_eqs())
    SymbolicUtils.simplify(input)
    assert_unchanged(input, snapshot)
end

@testset "simplify vs simplify! (noise)" begin
    base() = complete(_noise_eqs())
    nonmut = SymbolicUtils.simplify(base())
    input = base()
    mut = simplify!(input)
    assert_pair_identical(nonmut, mut)
end

# ---- modify_equations / modify_equations! -----------------------------------

@testset "modify_equations vs modify_equations! (deterministic)" begin
    @variables marker
    f = (op, rhs) -> rhs + marker
    nonmut = modify_equations(complete(_det_eqs()), f)
    input = complete(_det_eqs())
    mut = modify_equations!(input, f)
    assert_pair_identical(nonmut, mut)
end

@testset "modify_equations leaves its input untouched" begin
    @variables marker
    f = (op, rhs) -> rhs + marker
    input = complete(_det_eqs())
    snapshot = complete(_det_eqs())
    modify_equations(input, f)
    assert_unchanged(input, snapshot)
end

@testset "modify_equations vs modify_equations! (noise)" begin
    @variables marker
    f = (op, rhs) -> rhs + marker
    nonmut = modify_equations(complete(_noise_eqs()), f)
    input = complete(_noise_eqs())
    mut = modify_equations!(input, f)
    assert_pair_identical(nonmut, mut)
end

# ---- substitute / substitute! -------------------------------------------------

@testset "substitute vs substitute! (deterministic)" begin
    @variables ω U marker
    dict = Dict(ω => marker, U => 0)
    nonmut = SymbolicUtils.substitute(complete(_det_eqs()), dict)
    input = complete(_det_eqs())
    mut = substitute!(input, dict)
    @test mut === input
    assert_pair_identical(nonmut, mut)
end

@testset "substitute leaves its input untouched" begin
    @variables ω marker
    input = complete(_det_eqs())
    snapshot = complete(_det_eqs())
    SymbolicUtils.substitute(input, Dict(ω => marker))
    assert_unchanged(input, snapshot)
end

@testset "substitute vs substitute! (noise)" begin
    @variables κ marker
    dict = Dict(κ => marker)
    nonmut = SymbolicUtils.substitute(complete(_noise_eqs()), dict)
    @test nonmut isa NoiseMeanfieldEquations
    input = complete(_noise_eqs())
    mut = substitute!(input, dict)
    assert_pair_identical(nonmut, mut)
end
