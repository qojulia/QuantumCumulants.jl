using QuantumCumulants
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports

@testset "Quality gates" begin
    @testset "Aqua" begin
        Aqua.test_all(
            QuantumCumulants;
            # The Aqua persistent-task probe shells out to a precompile sandbox
            # for the full dependency stack (MTKBase, Symbolics, ...); on this
            # rewrite branch it routinely exceeds Aqua's 5-minute timeout, even
            # though the package itself has no persistent tasks.
            persistent_tasks = false,
            stale_deps = (; ignore = [:LaTeXStrings,
                                       :Latexify, :MacroTools, :OrderedCollections,
                                       :TermInterface]),
            piracies = (;
                # Base.complex glue methods are intentionally pirated to bridge
                # MTK's runtime code-gen for SQA-generated equations.
                # SecondQuantizedAlgebra.average(op, order; ...) is a deliberate
                # convenience overload that immediately applies cumulant_expansion;
                # the method takes only built-in types so Aqua flags it.
                treat_as_own = [Base.complex,
                                QuantumCumulants.SecondQuantizedAlgebra.average],
            ),
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
