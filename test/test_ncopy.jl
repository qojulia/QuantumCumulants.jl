using Qumulants
using SymbolicUtils

# Hilbert space
N = 10
hf = FockSpace(:field)
ha = NLevelSpace(:atom,(:g,:e);n=N)
# hN = ncopy(ha, N)
h_tot = hf⊗ha

# Operators
a = Destroy(h_tot,:a)
σ(i,j,k) = Transition(h_tot, :σ, i, j; index=k)

σ(k) = σ(:g,:e,k)

H = +((σ(i)'*σ(i) for i=1:N)...)

he = heisenberg(σ(1), H)


swap_idx(σ(1), 2)

swap_idx(σ(1)*σ(2), )
