using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test

const QC = QuantumCumulants

Base.@constprop :none function allocated_call(f, args::Vararg{Any, N}) where {N}
    b0 = Ref{Int64}(0)
    b1 = Ref{Int64}(0)
    Base.gc_bytes(b0)
    Base.@noinline f(args...)
    Base.gc_bytes(b1)
    return b1[] - b0[]
end

function numeric_coefficients(ir, pdict)
    pd = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in pdict)
    for coeff in ir.coeffs
        coeff isa Number && continue
        for var in Symbolics.get_variables(coeff)
            u = Symbolics.unwrap(var)
            if SymbolicUtils.issym(u) &&
                    Base.nameof(u) === :im &&
                    SymbolicUtils.symtype(u) === Number
                pd[u] = im
            end
        end
    end
    return ComplexF64[
        ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(coeff, pd))) for
            coeff in ir.coeffs
    ]
end

function reference_du(eqs, pdict, u)
    subs = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in pdict)
    for (i, state) in enumerate(eqs.states)
        subs[Symbolics.unwrap(state)] = u[i]
        adj = Symbolics.unwrap(average(adjoint(undo_average(state))))
        haskey(subs, adj) || (subs[adj] = conj(u[i]))
    end
    out = Vector{ComplexF64}(undef, length(eqs.states))
    for (i, eq) in enumerate(eqs.equations)
        rhs = Symbolics.unwrap(eq.rhs)
        for var in Symbolics.get_variables(rhs)
            uvar = Symbolics.unwrap(var)
            if SymbolicUtils.issym(uvar) &&
                    Base.nameof(uvar) === :im &&
                    SymbolicUtils.symtype(uvar) === Number
                subs[uvar] = im
            end
        end
        out[i] = ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(rhs, subs)))
    end
    return out
end

function pauli_fixture()
    n = 3
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:n]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hx γ
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(n - 1)) - hx * sum(σx(i) for i in 1:n)
    eqs = meanfield(
        [σz(i) for i in 1:n],
        H,
        [σm(i) for i in 1:n];
        rates = fill(γ, n),
        order = 2,
    )
    complete!(eqs)
    return eqs, Dict(J => 1.0, hx => 1.0, γ => 0.2)
end

@testset "serial MomentKernel agrees with symbolic substitution" begin
    eqs, ps = pauli_fixture()
    ir = QC._lower_moment_ir(eqs)
    kernel = QC.MomentKernel(ir)
    coeffs = numeric_coefficients(ir, ps)
    u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in eachindex(eqs.states)]
    du = similar(u)

    kernel(du, u, coeffs)
    ref = reference_du(eqs, ps, u)
    @test maximum(abs.(du .- ref) ./ max.(abs.(ref), 1.0e-12)) < 1.0e-12
end

@testset "serial MomentKernel handles folded conjugates" begin
    h = FockSpace(:kernel_cavity)
    a = Destroy(h, :a)
    @variables Δ Ω κ U
    H = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    complete!(eqs; get_adjoints = false)
    ps = Dict(Δ => -1.0, Ω => 1.3, κ => 1.0, U => 0.1)

    ir = QC._lower_moment_ir(eqs)
    kernel = QC.MomentKernel(ir)
    coeffs = numeric_coefficients(ir, ps)
    u = ComplexF64[0.3cos(2.1i) + 0.2im * sin(0.7i) for i in eachindex(eqs.states)]
    du = similar(u)

    kernel(du, u, coeffs)
    @test du ≈ reference_du(eqs, ps, u) rtol = 1.0e-12 atol = 1.0e-12
end

@testset "MomentKernel hot call is allocation-free after scratch initialization" begin
    eqs, ps = pauli_fixture()
    ir = QC._lower_moment_ir(eqs)
    kernel = QC.MomentKernel(ir)
    coeffs = numeric_coefficients(ir, ps)
    u = zeros(ComplexF64, length(eqs.states))
    du = similar(u)

    kernel(du, u, coeffs)
    allocated = allocated_call(kernel, du, u, coeffs)
    @test allocated == 0
end

@testset "MomentKernel is reentrant across concurrent tasks" begin
    eqs, ps = pauli_fixture()
    ir = QC._lower_moment_ir(eqs)
    kernel = QC.MomentKernel(ir)
    coeffs = numeric_coefficients(ir, ps)
    n = length(eqs.states)

    us = [
        ComplexF64[0.1cos(2.3i + 0.7j) + 0.05im * sin(1.1i - 0.3j) for i in 1:n] for
            j in 1:8
    ]
    refs = map(us) do u
        du = similar(u)
        kernel(du, u, coeffs)
        du
    end

    mismatches = Threads.Atomic{Int}(0)
    @sync for _ in 1:100, (u, ref) in zip(us, refs)
        Threads.@spawn begin
            du = similar(u)
            kernel(du, u, coeffs)
            du == ref || Threads.atomic_add!(mismatches, 1)
        end
    end
    @test mismatches[] == 0
end

@testset "task-local scratch does not retain discarded kernels" begin
    function weak_kernel_after_call()
        eqs, ps = pauli_fixture()
        ir = QC._lower_moment_ir(eqs)
        kernel = QC.MomentKernel(ir)
        weak_kernel = Base.WeakRef(kernel)
        coeffs = numeric_coefficients(ir, ps)
        u = zeros(ComplexF64, length(eqs.states))
        du = similar(u)
        kernel(du, u, coeffs)
        kernel = nothing
        GC.gc()
        return weak_kernel
    end
    @test weak_kernel_after_call().value === nothing
end
