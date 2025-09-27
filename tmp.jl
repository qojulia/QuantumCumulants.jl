using ModelingToolkit
using QuantumCumulants
using OrdinaryDiffEq

order = 2 #order of the cumulant expansion
@cnumbers N Δ g κ Γ R ν M

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)

h = hc ⊗ ha

k = Index(h, :k, N, ha)
l = Index(h, :l, N, ha)
m = Index(h, :m, N, ha)
n = Index(h, :n, N, ha)

@qnumbers a::Destroy(h)

σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
σn(i, j, k) = NumberedOperator(Transition(h, :σ, i, j), k)

# Define the Hamiltonian
H = -Δ*a'a + g*(Σ(a'*σ(1, 2, k), k) + Σ(a*σ(2, 1, k), k))

J = [a, σ(1, 2, l), σ(2, 1, l), σ(2, 2, l)]
rates = [κ, Γ, R, ν]
ops = [a'*a, σ(2, 2, m)]
eqs = meanfield(ops, H, J; rates = rates, order = order)

# System parameters
N_ = 2e5
Γ_ = 1.0 #Γ=1mHz
Δ_ = 2500Γ_ #Δ=2.5Hz
g_ = 1000Γ_ #g=1Hz
κ_ = 5e6*Γ_ #κ=5kHz
R_ = 1000Γ_ #R=1Hz
ν_ = 1000Γ_ #ν=1Hz

ps = [N, Δ, g, κ, Γ, R, ν]
p0 = [N_, Δ_, g_, κ_, Γ_, R_, ν_]

op1 = σ(1, 2, m)
op2 = σ(2, 1, n)

# Test system with no filter func to evaluate corr function

eqs_c_2 = complete(eqs);
# eqs_ev = evaluate(eqs_c_2; limits = (N=>3))

# u0_ev = zeros(ComplexF64, length(eqs_ev))

# @named sys_ev = ODESystem(eqs_ev)
# prob_ev = ODEProblem(sys_ev, u0_ev, (0.0, 0.1/50Γ_), ps .=> p0);
# sol_ev = solve(prob_ev, Tsit5(), maxiters = 1e7)

corr3 = CorrelationFunction(op1, op2, eqs_c_2)
corr3_ev = evaluate(corr3, 1, 2; limits = (N=>3));

@named csys3 = ODESystem(corr3_ev)
equations(csys3)[15]




# ps = []
# for eq ∈ corr3_ev.de.equations
#     ModelingToolkit.collect_vars!([], ps, eq.rhs, τ)
# end
# unique!(ps)
# @named csys3 = ODESystem(corr3_ev)

# c = corr3_ev
# τ = ModelingToolkit.get_iv(c.de)

# ps = QuantumCumulants.extract_parameters(corr3_ev.de) # [Δ,g,κ,Γ,ν,R]

# avg = average(c.op2_0) # ⟨σ212⟩
# avg2 = average(c.op2) # ⟨σ_0212⟩
# ps_ = [ps..., avg]
# de = substitute(c.de, Dict(avg2 => avg))

# ps_avg = filter(x->x isa Average, ps_) # ⟨σ212⟩
# ps_adj = map(QuantumCumulants._conj, ps_avg) # ⟨σ122⟩
# filter!(x->!(x ∈ Set(ps_avg)), ps_adj)
# ps_adj_hash = hash.(ps_adj)

# de_ = deepcopy(de)
# for i = 1:length(de.equations)
#     lhs = de_.equations[i].lhs
#     rhs = QuantumCumulants.substitute_conj(de_.equations[i].rhs, ps_adj, ps_adj_hash)
#     de_.equations[i] = Symbolics.Equation(lhs, rhs)
# end

# avg0 = average(c.op2_0) # ⟨σ212⟩
# avg0_par = QuantumCumulants._make_parameter(avg0) # var"⟨σ212⟩"
# de_ = substitute(de_, Dict(avg0 => avg0_par))

# eqs = ModelingToolkit.equations(de_)
# # eqs[15]
# for eq in eqs
#     println(eq)
# end
# sys = MTK.ODESystem(eqs, τ; kwargs...)
# return complete_sys ? complete(sys) : sys

# equations(csys3)
# unknowns(csys3)
# parameters(csys3)

val = ComplexF64[
    -1.1437146742985789e-7-6.3887991023161414e-24im,
    0.8635967567168613+0.0im,
    0.8635967567168613+0.0im,
    0.8635967567168613+0.0im,
    3.452543902275645e-7-0.00034541449685634095im,
    3.4525439022756444e-7-0.000345414496856341im,
    3.452543902275641e-7-0.0003454144968563411im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.00016965965013210335-2.727267523449695e-23im,
    0.0001696596501321033-3.242573084651477e-23im,
    0.00016965965013210333+4.61495050759439e-24im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.7457993582118818+0.0im,
    0.7457993582118818+3.805486687267425e-24im,
    0.7457993582118817+0.0im,
    0.0+0.0im,
    0.0+0.0im,
    0.0+0.0im,
]
p0_c3 = correlation_p0(corr3_ev, val, ps .=> p0)
u0_c3 = correlation_u0(corr3_ev, val)
# isequal.(first.(u0_c3),unknowns(csys3))

dict = merge(Dict(u0_c3), Dict(p0_c3))

# # "var\"⟨σ122⟩\"(τ)" ∈ string.(unknowns(csys3))
# # "var\"⟨σ122⟩\"(τ)" ∈ string.(first.(u0_c3))

# for eq in equations(csys3)
#     println(eq)
# end
# equations(csys3)[15]

prob_c3 = ODEProblem(csys3, dict, (0.0, 0.05); fully_determined = true)
# prob_c3.f.f.f_oop
sol_c3 = solve(prob_c3, Tsit5(); saveat = 0.001, maxiters = 1e8);

sys = csys3
op = dict
tspan = (0.0, 0.05)
callback = nothing
check_length = true
eval_expression = false, eval_module = @__MODULE__
check_compatibility = true, dvs = unknowns(sys)
ps = parameters(sys; initial_parameters = true)
iv = ModelingToolkit.has_iv(sys) ? ModelingToolkit.get_iv(sys) : nothing
eqs = equations(sys)

u0Type = pType = typeof(op)

op = ModelingToolkit.to_varmap(op, dvs)

ModelingToolkit.symbols_to_symbolics!(sys, op)
ModelingToolkit.check_inputmap_keys(sys, op)

defs = ModelingToolkit.add_toterms(
    ModelingToolkit.recursive_unwrap(defaults(sys));
    replace = ModelingToolkit.is_discrete_system(sys),
)
obs, eqs = ModelingToolkit.unhack_observed(observed(sys), eqs)


u0map = ModelingToolkit.anydict()
pmap = ModelingToolkit.anydict()
missing_unknowns, missing_pars =
    ModelingToolkit.build_operating_point!(sys, op, u0map, pmap, defs, dvs, ps)

u0_eltype = nothing
floatT = ModelingToolkit.calculate_float_type(op, u0Type)
u0_eltype = ModelingToolkit.something(u0_eltype, floatT)

u0_constructor = identity
p_constructor = identity
symbolic_u0 = false
u0_constructor =
    ModelingToolkit.get_u0_constructor(u0_constructor, u0Type, u0_eltype, symbolic_u0)
p_constructor = ModelingToolkit.get_p_constructor(p_constructor, pType, floatT)

guesses_ = ModelingToolkit.AnyDict()
t = zero(floatT)




has_u0_ics = false
op = copy(ModelingToolkit.anydict(op))
initialization_eqs = []
check_units = false
time_dependent_init = true
algebraic_only = false
isys = generate_initializesystem(
    sys;
    op,
    initialization_eqs,
    check_units,
    time_dependent_init,
    guesses = guesses_,
    algebraic_only,
)
simplify_system = true

uninit = setdiff(unknowns(sys), unknowns(isys), observables(isys))

# TODO: throw on uninitialized arrays
filter!(x -> !(x isa Symbolics.Arr), uninit)
