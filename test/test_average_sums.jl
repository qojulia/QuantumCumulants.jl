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

sum1 = SingleSum(σ(1,2,ind(:k)),ind(:k))
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
@test g_ isa SymbolicUtils.BasicSymbolic
@test gamma isa SymbolicUtils.BasicSymbolic{qc.DoubleNumberedVariable}

gamma_ = insert_index(gamma,ind(:j),2)
# @test gamma_ isa SymbolicUtils.BasicSymbolic{Parameter}
@test gamma_ isa SymbolicUtils.BasicSymbolic{Complex{Real}}

@test !isequal(gamma,gamma_)
@test !isequal(gamma,g_)
@test isequal(gamma2_,gamma)
@test !isequal(gamma,gamma2)

sumterm = σ(1,2,ind(:i))*σ(2,1,ind(:j))*σ(2,2,ind(:k))
sum_ = Σ(sumterm,ind(:i),[ind(:j),ind(:k)])
sum_A = average(sum_)

@test isequal(cumulant_expansion(sum_A,2),qc.IndexedAverageSum(cumulant_expansion(average(sumterm),2),ind(:i),[ind(:j),ind(:k)]))

inds = qc.get_indices(sumterm)
@test isequal([ind(:i),ind(:j),ind(:k)],inds)

pind = Index(h,:p,5,ha)
@test isequal(4*a,Σ(a,pind,[ind(:i)]))
@test isequal(average(Σ(σ(1,2,ind(:i)),ind(:i))),qc.IndexedAverageSum(average(σ(1,2,ind(:i))),ind(:i),[]))

avrgTerm = average(Σ(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i)))
@test avrgTerm isa SymbolicUtils.BasicSymbolic && operation(avrgTerm) === +
# ADsum1 = qc.IndexedAverageDoubleSum(avrgTerm,ind(:j),[ind(:i)])
# @test ADsum1 isa SymbolicUtils.BasicSymbolic && operation(ADsum1) === +
# @test SymbolicUtils.metadata(arguments(ADsum1)[it_2])[qc.IndexedAverageDoubleSum] isa qc.IndexedAverageDoubleSum
# @test SymbolicUtils.metadata(arguments(ADsum1)[2])[qc.IndexedAverageSum] isa qc.IndexedAverageSum

# executing it line by line gives a different order of the arguments?? arguments(ADsum1)
ADsum1 = simplify(qc.IndexedAverageDoubleSum(avrgTerm,ind(:j),[ind(:i)]))
it_1 = findfirst(x->typeof((x))==SymbolicUtils.BasicSymbolic{IndexedAverageSum}, arguments(ADsum1))
it_2 = findfirst(x->typeof((x))==SymbolicUtils.BasicSymbolic{IndexedAverageDoubleSum}, arguments(ADsum1))
@test ADsum1 isa SymbolicUtils.BasicSymbolic && operation(ADsum1) === +
@test SymbolicUtils.metadata(arguments(ADsum1)[it_2])[qc.IndexedAverageDoubleSum] isa qc.IndexedAverageDoubleSum
@test SymbolicUtils.metadata(arguments(ADsum1)[it_1])[qc.IndexedAverageSum] isa qc.IndexedAverageSum


@test isequal(qc.SpecialIndexedAverage(average(σ(1,2,ind(:i))),[(ind(:i),ind(:j))])+qc.SpecialIndexedAverage(average(σ(2,1,ind(:j))),[(ind(:i),ind(:j))]),
qc.SpecialIndexedAverage(average(σ(1,2,ind(:i))) + average(σ(2,1,ind(:j))),[(ind(:i),ind(:j))]))

@test qc.SpecialIndexedAverage(average(0),[(ind(:i),ind(:j))]) == 0
@test SymbolicUtils.metadata(qc.SpecialIndexedAverage(average(σ(2,1,ind(:i))),[(ind(:i),ind(:j))]))[qc.SpecialIndexedAverage] isa qc.SpecialIndexedAverage

@test qc.undo_average(arguments(ADsum1)[it_2]) isa qc.DoubleSum
@test isequal(simplify(Σ(Σ(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i)),ind(:j),[ind(:i)])),simplify(qc.undo_average(ADsum1)))

@test σ(1,2,ind(:i))*σ(2,1,ind(:j))*σn(2,2,3) isa qc.QMul
@test σn(2,2,3)*σ(1,2,ind(:i))*σ(2,1,ind(:j)) isa qc.QMul

# @test SymbolicUtils.iscall(sum_A.metadata) == false
@test qc.IndexedAverageSum(1) == 1

specAvrg = qc.SpecialIndexedAverage(average(σ(2,1,ind(:i))*σ(1,2,ind(:j))),[(ind(:i),ind(:j))])

@test isequal("(i≠1)",qc.SecondQuantizedAlgebra.writeNeqs([(ind(:i),1)]))
@test (isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[it_2]),SymbolicUtils.arguments(avrgTerm)[1]) || isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[it_2]),SymbolicUtils.arguments(avrgTerm)[2]) )
# SymbolicUtilsv1.4.0 argument order changed
@test isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[it_2])),SymbolicUtils.arguments(average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))
@test isequal(SymbolicUtils.arguments(specAvrg),SymbolicUtils.arguments(average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))

@test isequal(qc.insert_index(σ(1,2,ind(:j))*σn(1,2,2),ind(:j),1),qc.insert_index(qc.insert_index(σ(1,2,ind(:i))*σ(1,2,ind(:j)),ind(:i),2),ind(:j),1))

dict_ = qc.create_value_map(g(ind(:i)),2)
dict = Dict{SymbolicUtils.BasicSymbolic,ComplexF64}()
push!(dict,(qc.SingleNumberedVariable(:g,1) => 2))
push!(dict,(qc.SingleNumberedVariable(:g,2) => 2))
@test isequal(dict_,dict)

@test isequal(qc.getAvrgs(specAvrg),average(σ(2,1,ind(:i))*σ(1,2,ind(:j))))
@test (isequal(qc.getAvrgs(SymbolicUtils.arguments(avrgTerm)[1]),average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))) || isequal(qc.getAvrgs(SymbolicUtils.arguments(avrgTerm)[2]),average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))
# SymbolicUtilsv1.4.0 argument order changed

@test isequal(qc.IndexedAverageSum(g(ind(:i)),ind(:i),[]),average(Σ(g(ind(:i)),ind(:i),[])))
@test qc.IndexedAverageSum(g(ind(:i)),ind(:i),[]) isa SymbolicUtils.BasicSymbolic{qc.IndexedAverageSum}

@test ind(:i) ∈ qc.get_indices(g(ind(:i)))

@cnumbers N_
ind2(i) = Index(h,i,N_,ha)

N_n = 10
mappingDict = Dict{SymbolicUtils.BasicSymbolic,Int64}(N_ => N_n)
sum2_A = average(∑(σ(1,2,ind2(:i))*σ(2,1,ind2(:j)),ind2(:i)))
@test isequal(qc.eval_term(sum2_A;limits=mappingDict),evaluate(sum2_A;limits=(N_ => N_n)))

@cnumbers n_ m_
ind3(i) = Index(h,i,(n_*m_),ha)
map2 = Dict{SymbolicUtils.BasicSymbolic,Int64}(n_ => 2, m_ => 2)
sum3_A = average(∑(σ(2,1,ind3(:i))*σ(1,2,ind3(:j)),ind3(:i)))
sum3_B = qc.insert_index(sum3_A,ind3(:j),2)
eva = qc.eval_term(sum3_B;limits=map2)
@test eva isa SymbolicUtils.BasicSymbolic{<:CNumber} && operation(eva) === +
@test length(arguments(eva)) == 4

@test qc.SecondQuantizedAlgebra.containsIndexedOps(average(a*σ(2,1,ind(:i))*σ(1,2,ind(:j))))
@test !(qc.SecondQuantizedAlgebra.containsIndexedOps(average(a'*a)))

end
