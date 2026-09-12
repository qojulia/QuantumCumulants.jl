using QuantumCumulants
using Test
using JET

const JET_CALL_THUNKS = Pair{String, Function}[]
const JET_OPT_THUNKS = Pair{String, Function}[]

function _jet_meanfield_simple()
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    return meanfield([a], H, [a]; rates = [κ])
end

push!(JET_CALL_THUNKS, "meanfield(damped cavity)" => _jet_meanfield_simple)

function _jet_cumulant_expand()
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    return cumulant_expansion(average(a' * a' * a), 2)
end

function _jet_complete_jc()
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    return complete!(eqs)
end

push!(JET_CALL_THUNKS, "cumulant_expansion" => _jet_cumulant_expand)
push!(JET_CALL_THUNKS, "complete!(JC)" => _jet_complete_jc)

function _jet_to_system_cavity()
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    complete!(eqs)
    return System(eqs; name = :cav)
end

push!(JET_CALL_THUNKS, "System(damped cavity)" => _jet_to_system_cavity)

const _JET_MOMENT_IR = QuantumCumulants.MomentIR(
    Any[:u],
    Int32[0, 1],
    Int32[0, 1],
    Int32[1, 2],
    Int32[2],
    Int32[1],
    Any[1],
    Any[],
)
const _JET_MOMENT_KERNEL = QuantumCumulants.MomentKernel(_JET_MOMENT_IR)
const _JET_MOMENT_U = ComplexF64[0.25 + 0.1im]
const _JET_MOMENT_COEFFS = ComplexF64[2.0 - 0.5im]
const _JET_KERNEL_PLAN = QuantumCumulants.KernelParameterPlan(_JET_MOMENT_IR)
const _JET_KERNEL_PARAMETERS = QuantumCumulants.KernelParameters(
    ComplexF64[],
    _JET_MOMENT_COEFFS,
)
const _JET_KERNEL_RHS = QuantumCumulants.KernelRHS(_JET_MOMENT_KERNEL, _JET_KERNEL_PLAN)

function _jet_moment_kernel_rhs()
    du = similar(_JET_MOMENT_U)
    _JET_MOMENT_KERNEL(du, _JET_MOMENT_U, _JET_MOMENT_COEFFS)
    return du
end

function _jet_kernel_rhs()
    du = similar(_JET_MOMENT_U)
    _JET_KERNEL_RHS(du, _JET_MOMENT_U, _JET_KERNEL_PARAMETERS, 0.0)
    return du
end

push!(JET_OPT_THUNKS, "MomentKernel hot RHS" => _jet_moment_kernel_rhs)
push!(JET_OPT_THUNKS, "KernelRHS hot call" => _jet_kernel_rhs)

@testset "Type Stability (JET)" begin
    @static if isempty(VERSION.prerelease)
        @testset "report_package" begin
            result = JET.report_package(
                QuantumCumulants;
                target_modules = (QuantumCumulants,),
                ignore_missing_comparison = true,
            )
            @test isempty(JET.get_reports(result))
        end

        @testset "report_call: $(name)" for (name, thunk) in JET_CALL_THUNKS
            rep = JET.@report_call target_modules = (QuantumCumulants,) ignore_missing_comparison = true thunk()
            @test isempty(JET.get_reports(rep))
        end

        @testset "report_opt: $(name)" for (name, thunk) in JET_OPT_THUNKS
            rep = JET.@report_opt target_modules = (QuantumCumulants,) thunk()
            @test isempty(JET.get_reports(rep))
        end
    end
end
