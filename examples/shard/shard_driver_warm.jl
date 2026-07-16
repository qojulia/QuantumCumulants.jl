# Warm-driver comparison for the RGF gate: does removing invokelatest actually buy warm
# runtime, and does dynamic dispatch over heterogeneous RGF types eat the win (fixed by
# FunctionWrappers)? Warm timings only, so fresh-session discipline is not required.
#   ORDER = 1 | 2 | 3 (default 3), CHUNK (default 10)

using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns, equations, parameters, get_iv
using Symbolics
using OrdinaryDiffEqLowOrderRK
using QuantumOptics
using RuntimeGeneratedFunctions
using FunctionWrappers: FunctionWrapper
RuntimeGeneratedFunctions.init(@__MODULE__)

isdefined(Main, :ORDER) || (ORDER = 3)
isdefined(Main, :CHUNK) || (CHUNK = 10)
const N = 6

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = ORDER)
complete!(eqs)
sys = mtkcompile(System(eqs; name = :sys))
dvs = unknowns(sys); ps = parameters(sys); tiv = get_iv(sys)
rhss = [eq.rhs for eq in equations(sys)]
neq = length(rhss)

ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:N]...)
u0 = initial_values(eqs, ψ0)
pvec = Float64[Dict(J => 1.0, hx => 1.0, γ => 0.2)[x] for x in ps]
tspan = (0.0, 5.0)

ranges = [r[1]:r[end] for r in Iterators.partition(1:neq, CHUNK)]
fexprs = [
    Symbolics.build_function(rhss[r], dvs, ps, tiv; expression = Val{true}, cse = true)[2]
        for r in ranges
]
fs = [@RuntimeGeneratedFunction(fe) for fe in fexprs]

const VT = SubArray{ComplexF64, 1, Vector{ComplexF64}, Tuple{UnitRange{Int64}}, true}
const FW = FunctionWrapper{Nothing, Tuple{VT, Vector{ComplexF64}, Vector{Float64}, Float64}}
fws = FW[
    FW(
            let g = f
                (du, u, p, t) -> (g(du, u, p, t); nothing)
        end
        ) for f in fs
]

driver_dyn = let fs = fs, ranges = ranges
    (du, u, p, t) -> (
        for (f, r) in zip(fs, ranges)
            f(view(du, r), u, p, t)
        end; nothing
    )
end
driver_fw = let fws = fws, ranges = ranges
    (du, u, p, t) -> (
        @inbounds for i in eachindex(fws)
            fws[i](view(du, ranges[i]), u, p, t)
        end; nothing
    )
end
driver_il = let fs = fs, ranges = ranges
    (du, u, p, t) -> (
        for (f, r) in zip(fs, ranges)
            Base.invokelatest(f, view(du, r), u, p, t)
        end; nothing
    )
end

du = zeros(ComplexF64, neq)
sols = Dict{Symbol, Any}()
for (name, drv) in ((:dyn, driver_dyn), (:fw, driver_fw), (:il, driver_il))
    drv(du, u0, pvec, 0.0)                              # compile
    tcall = minimum(@elapsed(drv(du, u0, pvec, 0.0)) for _ in 1:200)
    prob = ODEProblem(ODEFunction(drv), u0, tspan, pvec)
    solve(prob, RK4(); saveat = 0.5)                    # compile solver path
    tsolve = minimum(@elapsed(sols[name] = solve(prob, RK4(); saveat = 0.5)) for _ in 1:3)
    println("DRIVER $name percall=$(round(tcall * 1.0e6, digits = 1))µs warmsolve=$(round(tsolve * 1000, digits = 1))ms")
end
d1 = maximum(maximum(abs, a - b) for (a, b) in zip(sols[:dyn].u, sols[:fw].u))
d2 = maximum(maximum(abs, a - b) for (a, b) in zip(sols[:dyn].u, sols[:il].u))
println("AGREEMENT dyn-fw=$(d1) dyn-il=$(d2)")
