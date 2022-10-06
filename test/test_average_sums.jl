using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics

const qc = QuantumCumulants

@testset "average_sums" begin

N = 2
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

ind(i) = Index(h,i,N,ha)

g(k) = IndexedVariable(:g,k)
Γij = DoubleIndexedVariable(:Γ,ind(:i),ind(:j),true)
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

@test(isequal(average(2*σ(1,2,ind(:k))),2*average(σ(1,2,ind(:k)))))
@test(isequal(average(g(ind(:k))*σ(2,2,ind(:k))),g(ind(:k))*average(σ(2,2,ind(:k)))))
@test(isequal(average(g(ind(:k))),g(ind(:k))))

sum1 = IndexedSingleSum(σ(1,2,ind(:k)),ind(:k))
σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
@test(isequal(evalTerm(average(sum1)),average(σn(1,2,1)) + average(σn(1,2,2))))
@test(isequal(σn(1,2,1)+σn(2,1,1),NumberedOperator(Transition(h,:σ,1,2)+Transition(h,:σ,2,1),1)))

#test insertIndex
@test(isequal(σn(2,2,1),insertIndex(σ(2,2,ind(:j)),ind(:j),1)))
@test(isequal(σ(1,2,ind(:j)),insertIndex(σ(1,2,ind(:j)),ind(:k),2)))
@test(isequal(1,insertIndex(1,ind(:k),1)))

sum2 = average(sum1*σ(1,2,ind(:l)))

@test(!isequal(σn(2,2,1),insertIndex(sum2,ind(:j),1)))


gamma = insertIndex(Γij,ind(:i),1)
@test insertIndex(g(ind(:j)),ind(:j),1) isa SymbolicUtils.Sym
@test gamma isa SymbolicUtils.Sym{Parameter,qc.numberedVariable}

@test insertIndex(gamma,ind(:j),2) isa SymbolicUtils.Sym

sumterm = σ(1,2,ind(:i))*σ(2,1,ind(:j))*σ(2,2,ind(:k))
sum_ = Σ(sumterm,ind(:i),[ind(:j),ind(:k)])
sum_A = average(sum_)

@test isequal(cumulant_expansion(sum_A,2),IndexedAverageSum(cumulant_expansion(average(sumterm),2),ind(:i),[ind(:j),ind(:k)]))

inds = qc.getIndices(sumterm)
@test isequal([ind(:i),ind(:j),ind(:k)],inds)



end

