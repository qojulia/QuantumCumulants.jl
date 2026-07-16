# PrecompileTools workload (issue #294, spec finding 12): one minimal end-to-end pass
# (N = 2 Pauli chain at order 2, meanfield -> complete! -> ODEProblem -> one RHS call)
# covering the completion machinery, the kernel lowering, and the generic kernel loops.
# No solver dependency; the sharded path is deliberately absent (RGF compilation cannot
# usefully precompile per-model generated code).
PrecompileTools.@compile_workload begin
    # `SQA.ProductSpace` instead of `⊗`: the tensor function is owned by
    # QuantumInterface, which the ExplicitImports quality gate will not import from
    h = SQA.ProductSpace((SQA.PauliSpace(:s1), SQA.PauliSpace(:s2)))
    s(k, i) = SQA.Pauli(h, :σ, k, i)
    Jp, hp, gp = Symbolics.@variables Jp hp gp
    Hp = -Jp * s(3, 1) * s(3, 2) - hp * (s(1, 1) + s(1, 2))
    eqs = meanfield(
        [s(3, 1), s(3, 2)], Hp, [(s(1, 1) - 1im * s(2, 1)) / 2];
        rates = [gp], order = 2,
    )
    complete!(eqs)
    prob = SciMLBase.ODEProblem(
        eqs, zeros(ComplexF64, length(eqs.states)), (0.0, 1.0),
        Dict(Jp => 1.0, hp => 1.0, gp => 0.1),
    )
    du = zeros(ComplexF64, length(eqs.states))
    prob.f(du, prob.u0, prob.p, 0.0)
end
