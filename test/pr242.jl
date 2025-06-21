using QuantumCumulants
using SymbolicUtils

@cnumbers g
hf = FockSpace(:cavity)
ha1 = NLevelSpace(:atom1, 2)
ha2 = NLevelSpace(:atom2, 2)
h = hf ⊗ ha1 ⊗ ha2
a = Destroy(h, :a)
s1(i, j) = Transition(h, :s1, i, j, 2)
s2(i, j) = Transition(h, :s2, i, j, 3)
H = g*(a' + a)*(s1(2, 1) + s1(1, 2) + s2(2, 1) + s2(1, 2))

rhs_ = commutator(im*H, a'a)
rhs_avg = average(rhs_)
x1 = SymbolicUtils.simplify(rhs_avg)

x2 = arguments(arguments(arguments(x1)[2])[3])
