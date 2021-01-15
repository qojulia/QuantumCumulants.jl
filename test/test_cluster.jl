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

# Collapse operators
J = [a;S(:g,:e);S(:e,:g)]
rates = [κ;[γ for i=1:M];[ν for i=1:M]]

# Derive equation for average photon number
he_n = average(heisenberg(a'*a,H,J;rates=rates), M)

id_aons = map(acts_on, S(:g,:e))
int_aons = [1]
he_scaled = scale(he_n, id_aons, int_aons, [N])
