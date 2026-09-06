using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test

const QC = QuantumCumulants

@testset "recursive moment polynomial lowering" begin
    h = PauliSpace(:s)
    sx = Pauli(h, :σ, 1)
    sy = Pauli(h, :σ, 2)
    sz = Pauli(h, :σ, 3)
    sm = (sx - 1im * sy) / 2
    @variables Ω γ a b
    eqs = meanfield([sz], Ω * sx, [sm]; rates = [γ], order = 1)
    complete!(eqs; get_adjoints = true)

    vars, idx = QC.statevars_resolved(eqs)
    @test length(vars) >= 2
    x, y = vars[1], vars[2]
    state_cache() = IdDict{Any, Bool}()

    # Nested sums/products and powers stay in the direct recursive path rather than requiring
    # whole-expression expansion through Symbolics.polynomial_coeffs.
    expr = Symbolics.unwrap(a * (x + b * y) * (x - y)^2)
    poly = QC._compile_moment_polynomial(expr, idx, state_cache())
    @test poly !== nothing

    generic, residual = Symbolics.polynomial_coeffs(expr, vars)
    @test QC._iszero_part(residual)
    generic_terms = Dict{Tuple, Any}(
        Tuple(QC.monomial_factors(mono, idx)) => coeff for (mono, coeff) in generic
    )
    @test Set(keys(poly)) == Set(keys(generic_terms))
    for mono in keys(poly)
        @test _is_zero(poly[mono] - generic_terms[mono])
    end

    # State-free denominators and arbitrary state-free symbolic subtrees remain coefficients.
    divided = Symbolics.unwrap((x + y)^2 / a)
    @test QC._compile_moment_polynomial(divided, idx, state_cache()) !== nothing
    constant_subtree = QC._compile_moment_polynomial(Symbolics.unwrap(sin(a)), idx, state_cache())
    @test length(constant_subtree) == 1
    @test haskey(constant_subtree, ())

    # Genuinely non-polynomial state dependence declines the fast path so the equation-local
    # generic fallback can preserve the established diagnostics.
    @test QC._compile_moment_polynomial(Symbolics.unwrap(exp(x)), idx, state_cache()) === nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap(x^(-1)), idx, state_cache()) === nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap(x^y), idx, state_cache()) === nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap((x + 1) / y), idx, state_cache()) === nothing

    # Sparse-polynomial helper identities and cancellation are part of the lowering contract.
    jx, jy = Int32(idx[x]), Int32(idx[y])
    mx, my = (jx,), (jy,)
    @test QC._merge_factor_tuples((), mx) == mx
    @test QC._merge_factor_tuples(mx, ()) == mx
    @test QC._merge_factor_tuples(mx, my) == Tuple(sort(Int32[jx, jy]))

    empty_poly = Dict{Tuple, Any}()
    px = Dict{Tuple, Any}(mx => 3)
    c2 = Dict{Tuple, Any}(() => 2)
    @test isempty(QC._poly_mul(empty_poly, px))
    @test isempty(QC._poly_mul(px, empty_poly))
    @test QC._poly_mul(c2, px)[mx] == 6
    @test QC._poly_mul(px, c2)[mx] == 6
    @test QC._poly_scale(px, 1) === px
    @test isempty(QC._poly_scale(px, 0))
    @test QC._poly_pow(px, 0) == Dict{Tuple, Any}(() => 1)
    @test QC._poly_pow(px, 1) === px

    cancelled = Dict{Tuple, Any}(mx => 1)
    QC._poly_add_term!(cancelled, mx, -1)
    @test isempty(cancelled)
    @test QC._constant_poly_coeff(Dict{Tuple, Any}(mx => 1)) === nothing

    # The generic monomial decoder remains the equation-local compatibility reference.
    @test isempty(QC.monomial_factors(1, idx))
    @test QC.monomial_factors(Symbolics.unwrap(x * y), idx) == sort(Int32[jx, jy])
    @test QC.monomial_factors(Symbolics.unwrap(x^2), idx) == Int32[jx, jx]
end
