using Qumulants

M = 2

# Define parameters
@parameters Δ g γ κ ν N

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(Symbol(:atoms),(:g,:e))
hc = ClusterSpace(ha,N,M)
h = hf⊗hc

a = Destroy(h,:a)
S(i,j) = Transition(h,:σ,i, j)

H = Δ*a'*a + g*(a'*sum(S(:g,:e)) + a*sum(S(:e,:g)))
H_ = simplify_operators(H)

# Collapse operators
J = [a;S(:g,:e);S(:e,:g)]
rates = [κ;[γ for i=1:M];[ν for i=1:M]]

# Derive equation for average photon number
he_n = heisenberg(a'*a,H,J;rates=rates)

he_avg = average(he_n, M)

phase_invariant(x) = iszero(phase(x))
phase(avg::Average) = phase(avg.operator)
phase(op::Destroy) = -1
phase(op::Create) = 1
function phase(t::Transition)
    lvls = Qumulants.levels(t.hilbert, t.aon)
    i = findfirst(isequal(t.i), lvls)
    j = findfirst(isequal(t.j), lvls)
    if i < j
        -1
    elseif i > j
        1
    else
        0
    end
end
phase(op::OperatorTerm{<:typeof(*)}) = sum(phase(arg) for arg in op.arguments)

he_nophase = complete(he_avg;order=M,filter_func=phase_invariant)
