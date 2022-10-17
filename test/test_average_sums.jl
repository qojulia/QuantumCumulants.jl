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
Γij = DoubleIndexedVariable(:Γ,ind(:i),ind(:j))
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
Ω(i,j) = IndexedVariable(:Ω,i,j;identical=false)

@test Ω(ind(:i),ind(:i)) == 0

a = Destroy(h,:a)

@test(isequal(average(2*σ(1,2,ind(:k))),2*average(σ(1,2,ind(:k)))))
@test(isequal(average(g(ind(:k))*σ(2,2,ind(:k))),g(ind(:k))*average(σ(2,2,ind(:k)))))
@test(isequal(average(g(ind(:k))),g(ind(:k))))

sum1 = IndexedSingleSum(σ(1,2,ind(:k)),ind(:k))
σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
@test(isequal(eval_term(average(sum1)),average(σn(1,2,1)) + average(σn(1,2,2))))
@test(isequal(σn(1,2,1)+σn(2,1,1),NumberedOperator(Transition(h,:σ,1,2)+Transition(h,:σ,2,1),1)))

#test insert_index
@test(isequal(σn(2,2,1),insert_index(σ(2,2,ind(:j)),ind(:j),1)))
@test(isequal(σ(1,2,ind(:j)),insert_index(σ(1,2,ind(:j)),ind(:k),2)))
@test(isequal(1,insert_index(1,ind(:k),1)))

sum2 = average(sum1*σ(1,2,ind(:l)))

@test(!isequal(σn(2,2,1),insert_index(sum2,ind(:j),1)))


gamma = insert_index(Γij,ind(:i),1)
gamma2 = insert_index(Γij,ind(:j),2)
gamma2_ = insert_index(Γij,ind(:i),1)
g_ = insert_index(g(ind(:j)),ind(:j),1)
@test g_ isa SymbolicUtils.Sym
@test gamma isa SymbolicUtils.Sym{Parameter,qc.numberedVariable}

gamma_ = insert_index(gamma,ind(:j),2)
@test gamma_ isa SymbolicUtils.Sym

@test !isequal(gamma,gamma_)
@test !isequal(gamma,g_)
@test isequal(gamma2_,gamma)
@test !isequal(gamma,gamma2)

sumterm = σ(1,2,ind(:i))*σ(2,1,ind(:j))*σ(2,2,ind(:k))
sum_ = Σ(sumterm,ind(:i),[ind(:j),ind(:k)])
sum_A = average(sum_)

@test isequal(cumulant_expansion(sum_A,2),qc.IndexedAverageSum(cumulant_expansion(average(sumterm),2),ind(:i),[ind(:j),ind(:k)]))

inds = qc.getIndices(sumterm)
@test isequal([ind(:i),ind(:j),ind(:k)],inds)

pind = Index(h,:p,5,ha)
@test isequal(4*a,Σ(a,pind,[ind(:i)]))
@test isequal(average(Σ(σ(1,2,ind(:i)),ind(:i))),qc.IndexedAverageSum(average(σ(1,2,ind(:i))),ind(:i),[]))

@test isequal(0,qc.IndexedAverageSum(average(σ(2,1,ind(:i))*σ(2,1,ind(:i))),ind(:i),[]))

avrgTerm = average(Σ(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i)))
@test avrgTerm isa SymbolicUtils.Add
ADsum1 = qc.IndexedAverageDoubleSum(avrgTerm,ind(:j),[ind(:i)])
@test ADsum1 isa SymbolicUtils.Add
@test arguments(ADsum1)[1].metadata isa qc.IndexedAverageDoubleSum  
@test arguments(ADsum1)[2].metadata isa qc.IndexedAverageSum  

@test qc.NumberedIndexedAverage(average(σ(1,2,ind(:i))),[]) isa Average
@test qc.NumberedIndexedAverage(average(σ(1,2,ind(:j))),[(ind(:j),3)]).metadata isa qc.NumberedIndexedAverage
@test qc.NumberedIndexedAverage(average(0),[(ind(:j),3)]) == 0

@test isequal(qc.SpecialIndexedAverage(average(σ(1,2,ind(:i))),[(ind(:i),ind(:j))])+qc.SpecialIndexedAverage(average(σ(2,1,ind(:j))),[(ind(:i),ind(:j))]),
qc.SpecialIndexedAverage(average(σ(1,2,ind(:i))) + average(σ(2,1,ind(:j))),[(ind(:i),ind(:j))]))

@test qc.SpecialIndexedAverage(average(0),[(ind(:i),ind(:j))]) == 0
@test qc.SpecialIndexedAverage(average(σ(2,1,ind(:i))),[(ind(:i),ind(:j))]).metadata isa qc.SpecialIndexedAverage

@test qc.undo_average(arguments(ADsum1)[1]) isa qc.IndexedDoubleSum
@test isequal(Σ(Σ(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i)),ind(:j),[ind(:i)]),qc.undo_average(ADsum1))

@test σ(1,2,ind(:i))*σ(2,1,ind(:j))*σn(2,2,3) isa qc.QMul
@test σn(2,2,3)*σ(1,2,ind(:i))*σ(2,1,ind(:j)) isa qc.QMul

@test SymbolicUtils.istree(sum_A.metadata) == false
@test qc.IndexedAverageSum(1) == 1

specAvrg = qc.SpecialIndexedAverage(average(σ(2,1,ind(:i))*σ(1,2,ind(:j))),[(ind(:i),ind(:j))])

@test isequal("(i≠1)",qc.writeNeqs([(ind(:i),1)]))
@test isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[1]),SymbolicUtils.arguments(avrgTerm)[1])
@test isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[1])),SymbolicUtils.arguments(average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))
@test isequal(SymbolicUtils.arguments(specAvrg),SymbolicUtils.arguments(average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))

@test isequal(σ(1,2,ind(:j))*σn(1,2,2),qc.insert_index(σ(1,2,ind(:i))*σ(1,2,ind(:j)),ind(:i),2))

vec1 = qc.appendEach([1,2,3],[3,4,5])
vec2 = [[1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5], [3, 3], [3, 4], [3, 5]]
@test isequal(vec1,vec2)

dict_ = qc.create_value_map(g(ind(:i)),2)
dict = Dict{SymbolicUtils.Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
push!(dict,(qc.SingleNumberedVariable(:g,1) => 2))
push!(dict,(qc.SingleNumberedVariable(:g,2) => 2))
@test isequal(dict_,dict)

@test isequal(qc.getAvrgs(specAvrg),[average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))])

@test isequal(qc.getAvrgs(SymbolicUtils.arguments(avrgTerm)[1]),[average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))])

@test isequal(ind(:i).rangeN * 5, qc.IndexedAverageSum(5,ind(:i),[]))
@test isequal(qc.IndexedAverageSum(g(ind(:i)),ind(:i),[]),average(Σ(g(ind(:i)),ind(:i)),[]))
@test qc.IndexedAverageSum(g(ind(:i)),ind(:i),[]) isa SymbolicUtils.Sym{Parameter,qc.IndexedAverageSum}

@test ind(:i) ∈ qc.getIndices(g(ind(:i)))

end

