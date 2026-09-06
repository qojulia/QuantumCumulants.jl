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

    # Nested sums/products and powers stay in the direct recursive path rather than requiring
    # whole-expression expansion through Symbolics.polynomial_coeffs.
    expr = Symbolics.unwrap(a * (x + b * y) * (x - y)^2)
    poly = QC._compile_moment_polynomial(expr, idx, IdDict{Any, Bool}())
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

    # State-free denominators are coefficients; genuinely non-polynomial state dependence
    # declines the fast path so the equation-local generic fallback can preserve diagnostics.
    divided = Symbolics.unwrap((x + y)^2 / a)
    @test QC._compile_moment_polynomial(divided, idx, IdDict{Any, Bool}()) !== nothing
    @test QC._compile_moment_polynomial(
        Symbolics.unwrap(exp(x)), idx, IdDict{Any, Bool}()
    ) === nothing
    @test QC._compile_moment_polynomial(
        Symbolics.unwrap(x^(-1)), idx, IdDict{Any, Bool}()
    ) === nothing
end
