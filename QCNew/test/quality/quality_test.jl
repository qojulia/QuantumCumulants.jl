using QCNew
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports

@testset "Quality gates" begin
    @testset "Aqua" begin
        Aqua.test_all(
            QCNew;
            # The persistent-task probe precompiles the full dependency stack
            # (MTKBase, Symbolics, ...) and exceeds Aqua's 5-minute timeout; the
            # package itself has no persistent tasks.
            persistent_tasks = false,
        )
    end

    @testset "ExplicitImports" begin
        @test check_no_implicit_imports(QCNew) === nothing
        @test check_all_explicit_imports_via_owners(QCNew) === nothing
        @test check_no_stale_explicit_imports(QCNew) === nothing
        @test check_all_qualified_accesses_via_owners(
            QCNew; skip = (Base => Core,),
        ) === nothing
        @test check_no_self_qualified_accesses(QCNew) === nothing
    end

    @testset "CheckConcreteStructs" begin
        for name in names(QCNew; all = true)
            isdefined(QCNew, name) || continue
            T = getfield(QCNew, name)
            T isa Type || continue
            isabstracttype(T) && continue
            T isa UnionAll && continue
            isstructtype(T) || continue
            parentmodule(T) === QCNew || continue
            @testset "$name" begin
                @test all_concrete(T; verbose = false)
            end
        end
    end
end
