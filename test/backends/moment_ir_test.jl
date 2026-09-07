using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test

const QC = QuantumCumulants

function ir_factor_tuples(ir)
    factors = Vector{Tuple}(undef, length(ir.parent))
    factors[1] = ()
    for m in 2:length(ir.parent)
        p = Int(ir.parent[m])
        factors[m] = (factors[p]..., ir.leaf[m])
    end
    return factors
end

function oracle_monomial_factors(mono, idx)
    fs = Int32[]
    addfactor(f) = if SymbolicUtils.iscall(f) && SymbolicUtils.operation(f) === (^)
        base, exponent = SymbolicUtils.arguments(f)
        n = Int(SymbolicUtils.unwrap_const(exponent))
        append!(fs, fill(idx[base], n))
    else
        push!(fs, idx[f])
    end

    if mono isa Number || SymbolicUtils.isconst(mono)
        return ()
    elseif SymbolicUtils.iscall(mono) && SymbolicUtils.operation(mono) === (*)
        foreach(addfactor, SymbolicUtils.arguments(mono))
    else
        addfactor(mono)
    end
    sort!(fs)
    return Tuple(fs)
end

function oracle_equation_terms(drift, vars, idx)
    dict, residual = Symbolics.polynomial_coeffs(drift, vars)
    @test _is_zero(residual)
    out = Dict{Tuple, Any}()
    for (mono, coeff) in dict
        factors = oracle_monomial_factors(mono, idx)
        out[factors] = haskey(out, factors) ? out[factors] + coeff : coeff
    end
    return out
end

function ir_equation_terms(ir, row)
    factors = ir_factor_tuples(ir)
    out = Dict{Tuple, Any}()
    for k in ir.rowptr[row]:(ir.rowptr[row + 1] - 1)
        mono = factors[ir.monomial[k]]
        coeff = ir.coeffs[ir.coeff_id[k]]
        out[mono] = haskey(out, mono) ? out[mono] + coeff : coeff
    end
    return out
end

function assert_ir_matches_symbolics(eqs)
    ir = QC._lower_moment_ir(eqs)
    idx = QC._moment_state_indices(eqs)
    vars = collect(keys(idx))
    @test length(ir.states) == length(eqs.states)
    @test length(ir.rowptr) == length(eqs.states) + 1
    @test ir.rowptr[1] == 1
    @test ir.rowptr[end] == length(ir.monomial) + 1

    for i in eachindex(eqs.equations)
        got = ir_equation_terms(ir, i)
        ref = oracle_equation_terms(Symbolics.unwrap(eqs.equations[i].rhs), vars, idx)
        @test Set(keys(got)) == Set(keys(ref))
        for mono in keys(ref)
            @test _is_zero(got[mono] - ref[mono])
        end
    end
    return ir
end

@testset "MomentIR polynomial oracle — Pauli hierarchy" begin
    h = PauliSpace(:s)
    sx = Pauli(h, :σ, 1)
    sy = Pauli(h, :σ, 2)
    sz = Pauli(h, :σ, 3)
    sm = (sx - 1im * sy) / 2
    @variables Ω γ
    eqs = meanfield([sz], Ω * sx, [sm]; rates = [γ], order = 2)
    complete!(eqs; get_adjoints = true)

    ir = assert_ir_matches_symbolics(eqs)
    ir2 = QC._lower_moment_ir(eqs)
    @test isequal(ir.states, ir2.states)
    @test ir.parent == ir2.parent
    @test ir.leaf == ir2.leaf
    @test ir.rowptr == ir2.rowptr
    @test ir.monomial == ir2.monomial
    @test ir.coeff_id == ir2.coeff_id
    @test isequal(ir.coeffs, ir2.coeffs)
    @test isequal(ir.params, ir2.params)
end

@testset "MomentIR polynomial oracle — folded nonlinear Kerr hierarchy" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables Δ Ω κ U
    H = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    complete!(eqs; get_adjoints = false)

    ir = assert_ir_matches_symbolics(eqs)
    @test any(x -> x < 0, ir.leaf)
end

@testset "recursive grammar and explicit capability boundary" begin
    h = PauliSpace(:grammar)
    sx = Pauli(h, :σ, 1)
    sy = Pauli(h, :σ, 2)
    sz = Pauli(h, :σ, 3)
    sm = (sx - 1im * sy) / 2
    @variables Ω γ a b
    eqs = meanfield([sz], Ω * sx, [sm]; rates = [γ], order = 1)
    complete!(eqs; get_adjoints = true)

    idx = QC._moment_state_indices(eqs)
    vars = collect(keys(idx))
    x, y = vars[1], vars[2]
    cache() = IdDict{Any, Bool}()

    nested = Symbolics.unwrap(a * (x + b * y) * (x - y)^2)
    @test QC._compile_moment_polynomial(nested, idx, cache()) !== nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap((x + y)^2 / a), idx, cache()) !== nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap(sin(a)), idx, cache()) !== nothing

    @test QC._compile_moment_polynomial(Symbolics.unwrap(exp(x)), idx, cache()) === nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap(x^(-1)), idx, cache()) === nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap(x^y), idx, cache()) === nothing
    @test QC._compile_moment_polynomial(Symbolics.unwrap((x + 1) / y), idx, cache()) === nothing

    nonpoly = modify_equations(eqs, (op, drift) -> drift + exp(average(op)))
    @test_throws QC.NonPolynomialDriftError QC._lower_moment_ir(nonpoly)

    timed = modify_equations(eqs, (op, drift) -> cos(eqs.iv) * drift)
    @test_throws QC.TimeDependentCoefficientError QC._lower_moment_ir(timed)
end

@testset "user parameter named im is distinct from the algebraic imaginary constant" begin
    h = PauliSpace(:im_parameter)
    sx = Pauli(h, :σ, 1)
    sy = Pauli(h, :σ, 2)
    sz = Pauli(h, :σ, 3)
    sm = (sx - 1im * sy) / 2
    @variables Ω γ
    eqs = meanfield([sz], Ω * sx, [sm]; rates = [γ], order = 1)
    complete!(eqs; get_adjoints = true)

    imvar = Symbolics.variable(:im)
    modified = modify_equations(eqs, (op, drift) -> imvar * drift)
    ir = QC._lower_moment_ir(modified)
    @test any(p -> isequal(p, Symbolics.unwrap(imvar)), ir.params)
end
