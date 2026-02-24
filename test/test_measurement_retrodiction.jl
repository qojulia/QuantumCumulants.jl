using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq
using StochasticDiffEq
using DiffEqNoiseProcess
using Random
using Test

# @testset "test_measurement_retrodiction" begin

    h = PhaseSpace(:motion)
    @qnumbers x::Position(h, 1)
    p = Momentum(h, :p, 1)
    @rnumbers Ω Γ η s
    m = 1
    @syms t::Real
    @register_symbolic Ydot(t)
    @register_symbolic Ydot_rev(t)
    a = (x + 1im*p)*s

    H = p^2/(2m) + 0.5m*Ω^2*x^2
    J = [a]
    R = [Γ]
    eff = [η]
    ops = [x, p, x*x, x*p, p*p]
    eqs = meanfield(ops, H, J; rates = R, efficiencies = eff, order = 2)

    Ω_ = 1
    Γ_ = 1/6
    η_ = 0.5
    Tend = 0.1/Γ_
    dt = Tend/5e3 # time step
    T_saveat = [0:dt:Tend;]
    ps = [Ω, Γ, η, s]
    pn = [Ω_, Γ_, η_, 1/√2]
    ## initial state
    x0 = 5.0;
    p0 = 0.0
    u0 = [x0, p0, 5+x0*x0, im/2+x0*p0, 5+p0*p0] # σx = σpp = 5

    Random.seed!(1)
    @named sys_fw = SDESystem(eqs)
    dict_fw = merge(Dict(unknowns(sys_fw) .=> u0), Dict(ps .=> pn))
    noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0; save_everystep = true)
    prob_fw = SDEProblem(sys_fw, dict_fw, (0, Tend); noise = noise)

    sol_fw = solve(prob_fw, EM(); dt = dt, saveat = T_saveat)

    t_W = noise.t # noise time steps
    W_W = noise.W # noise values
    N = length(t_W)-1
    dY_W = zeros(N)
    sol_x = real(sol_fw[x])
    for k = 1:N
        dW = W_W[k+1] - W_W[k]
        cd_p_c = sol_x[k]*√(2)*√(Γ_)

        dY_W[k] = dW + √η_*cd_p_c*dt
    end
    dYdt_data = [dY_W; dY_W[end]] / dt

    # dYdt(t) = dYdt_data[Int(floor(t/dt)) + 1] # TODO: why does this give a small difference?? Maybe the times do no perfectly match in the solver?
    dYdt(t) = dYdt_data[Int(round(t/dt))+1]
    Ydot(t) = dYdt(t)

    dYdt_data_rev = [reverse(dY_W); dY_W[1]] / dt
    dYdt_rev(t) = dYdt_data_rev[Int(round(t/dt))+1]
    Ydot_rev(t) = dYdt_rev(t)

    # modify equations
    eqs_det = meanfield(ops, H, J; rates = R, order = 2)
    function f_measure(lhs, rhs)
        term_ = √(η*Γ)*(lhs*a + a'*lhs - average(a+a')*lhs)*(Ydot(t) - √(η*Γ)*average(a+a'))
        term = cumulant_expansion(average(term_), 2)
        return rhs + term
    end
    eqs_kal = modify_equations(eqs_det, f_measure)
    @test isequal(eqs_kal[1].lhs, eqs_det[1].lhs)
    eqs_kal_1_rhs = eqs_det[1].rhs -s*(Ydot(t) - 2s*sqrt(Γ*η)*average(x))*sqrt(Γ*η) + 2s*(Ydot(t) - 2s*sqrt(Γ*η)*average(x))*sqrt(Γ*η)*average(x*x) - 2s*(Ydot(t) - 2s*sqrt(Γ*η)*average(x))*sqrt(Γ*η)*(average(x)^2)
    @test isequal(eqs_kal[1].rhs - eqs_kal_1_rhs, 0)
    eqs_kal_test = deepcopy(eqs_det)
    modify_equations!(eqs_kal_test, f_measure)
    @test isequal(eqs_kal_test[1], eqs_kal[1])
    @test isequal(eqs_kal_test[4], eqs_kal[4])

    @named sys_kal = System(eqs_kal)
    prob_kal = ODEProblem(sys_kal, dict_fw, (0, Tend))

    sol_kal = solve(prob_kal, Euler(); dt = dt, saveat = T_saveat)
    # it = 100#8 # step from 7 to 8 changes drastically if floor() is used!! # TODO
    @test sum(abs.(sol_fw.u[end] .- sol_kal.u[end])) < 1e-8

    ### back propagation
    eqs_back_kal = meanfield(ops, -H, J; rates = R, order = 2)
    eqs_back_kal_c = complete(eqs_back_kal)
    ## adjust recycling term for the back propagation
    function f_back_lind(lhs, rhs)
        term_1 = Γ*a*lhs*a' - Γ*a'lhs*a # adapt recycling term
        term_2 = -Γ*(average(a*a') - average(a'a))*lhs
        term = cumulant_expansion(average(term_1 + term_2), 2)
        return rhs + term
    end
    eqs_back_kal_c1 = modify_equations(eqs_back_kal_c, f_back_lind)
    ## include measurement backaction 
    function f_back_kal(lhs, rhs)
        term_ =
            √(η*Γ)*(lhs*a' + a*lhs - average(a+a')*lhs)*(Ydot_rev(t) - √(η*Γ)*average(a+a'))
        term = cumulant_expansion(average(term_), 2)
        return rhs + term
    end
    eqs_back_kal_c2 = modify_equations(eqs_back_kal_c1, f_back_kal)
    # deterministic smoothing equations
    # Ydot_rev(t) = dYdt(Tend-t)
    @named sys_back_kal = System(eqs_back_kal_c2)
    u0_bw = [0.0+0im, 0, 100, 0, 100]
    dict_bw = merge(Dict(unknowns(sys_back_kal) .=> u0_bw), Dict(ps .=> pn))
    prob_bw_kal = ODEProblem(sys_back_kal, dict_bw, (0, Tend))
    sol_bw = solve(prob_bw_kal, Euler(); dt = dt, saveat = T_saveat)

    # TODO: stochastic version - compare and test functions; 

    ##########################
    ### stochastic verison ### (comparison)
    ##########################

    # equations for the measuremt record dY as noise input for the solver (instead of dW)
    eqs_c_Y = translate_W_to_Y(eqs)

    # solve forward propagation again 
    Y_W_fw = cumsum([0; dY_W]) 
    noise_Y_fw = NoiseGrid(t_W, Y_W_fw)

    @named sys_fw_Y = SDESystem(eqs_c_Y)
    prob_fw_Y = SDEProblem(sys_fw_Y, dict_fw, (0, Tend); noise=noise_Y_fw)
    sol_fw_Y = solve(prob_fw_Y, EM(); dt=dt, saveat=T_saveat)
    @test sum(abs.(sol_fw.u[end] .- sol_fw_Y.u[end])) < 1e-8

    # backward equations
    eqs_back = meanfield_backward(ops[1:end-1], H, J; rates=R, efficiencies=eff, order=2)
    eqs_back_c = complete(eqs_back)
    @test isequal(eqs_back_c.states, eqs_c_Y.states)

    Y_W_bw = cumsum([0; reverse(dY_W)]) 
    noise_Y_bw = NoiseGrid(t_W, Y_W_bw)

    # noise_Y_bw = NoiseGrid(t_W, reverse(-Y_W))
    @named sys_fw_Y_bw = SDESystem(eqs_back_c)
    dict_bw = merge(Dict(unknowns(sys_fw_Y_bw) .=> u0_bw), Dict(ps.=>pn))
    prob_fw_Y_bw = SDEProblem(sys_fw_Y_bw, dict_bw, (0, Tend); noise=noise_Y_bw)
    sol_bw_Y = solve(prob_fw_Y_bw, EM(); dt=dt, saveat=T_saveat)

    @test sum(abs.(sol_bw.u[end] .- sol_bw_Y.u[end])) < 1e-8
# end
