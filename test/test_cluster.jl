using QuantumCumulants
using Test

@testset "cluster" begin

order = 4
N_c = 3 #number of clusters
N = [Parameter(Symbol(:N_, i)) for i=1:N_c]

hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atoms, j),3) for j=1:N_c]
ha_c = [ClusterSpace(ha[j],N[j],order) for j=1:N_c]
h = ⊗(hf, ha_c...)
# Define the fundamental operators
a = Destroy(h,:a,1)
S(i,j,c) = Transition(h,Symbol(:σ, c),i, j, 1+c) #c=cluster

@test QuantumCumulants.has_cluster(h)
@test !(QuantumCumulants.has_cluster(hf))

@test S(2,2,1) ≠ S(2,2,2) ≠ S(2,2,3)
@test isequal(S(2,2,2)[1]*S(2,2,1)[2], S(2,2,1)[2]*S(2,2,2)[1])
@test length(S(2,2,1)) == length(S(2,2,2)) == length(S(2,2,3)) == order
@test acts_on(S(2,2,1)[1]) == QuantumCumulants.ClusterAon(2,1)
@test acts_on(S(2,2,2)[2]) == QuantumCumulants.ClusterAon(3,2)
@test acts_on(S(2,2,1)[1]) < acts_on(S(2,2,1)[2]) < acts_on(S(2,2,2)[1])

end #testset
