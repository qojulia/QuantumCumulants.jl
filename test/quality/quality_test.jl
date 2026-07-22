using QuantumCumulants
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports

@testset "Quality gates" begin
    @testset "Aqua" begin
        Aqua.test_all(
            QuantumCumulants;
            # The persistent-task probe precompiles the full dependency stack
            # (MTKBase, Symbolics, ...) and exceeds Aqua's 5-minute timeout; the
            # package itself has no persistent tasks.
            persistent_tasks = false,
        )
    end

    @testset "ExplicitImports" begin
        @test check_no_implicit_imports(QuantumCumulants) === nothing
        @test check_all_explicit_imports_via_owners(QuantumCumulants) === nothing
        @test check_no_stale_explicit_imports(QuantumCumulants) === nothing
        @test check_all_qualified_accesses_via_owners(
            QuantumCumulants; skip = (Base => Core,),
        ) === nothing
        @test check_no_self_qualified_accesses(QuantumCumulants) === nothing

        @test check_all_explicit_imports_are_public(QuantumCumulants) === nothing
        @test check_all_qualified_accesses_are_public(
            QuantumCumulants;
            ignore = (
                :FnType, # SU.jl
                :isconst, # SU.jl
                :IM, # Symbolics.jl
                :RefValue, # Base
                :tobrownian, # MTK.jl
                :toparam # MTK.jl
                ),
        ) === nothing
    end

    @testset "CheckConcreteStructs" begin
        for name in names(QuantumCumulants; all = true)
            isdefined(QuantumCumulants, name) || continue
            T = getfield(QuantumCumulants, name)
            T isa Type || continue
            isabstracttype(T) && continue
            T isa UnionAll && continue
            isstructtype(T) || continue
            parentmodule(T) === QuantumCumulants || continue
            @testset "$name" begin
                @test all_concrete(T; verbose = false)
            end
        end
    end
end
