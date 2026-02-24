using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq 
using StochasticDiffEq
using Random 
using Test

@testset "test_measurement_retrodiction" begin

    h = PhaseSpace(:motion)
    @qnumbers x::Position(h,1)
    p = Momentum(h,:p,1)
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
    ops = [x,p,x*x,x*p,p*p] 
    eqs = meanfield(ops, H, J; rates=R, efficiencies=eff, order=2)

    Ω_ = 1
    Γ_ = 1/6
    η_ = 0.5
    Tend = 0.1/Γ_
    dt = Tend/5e3 # time step
    T_saveat = [0:dt:Tend;] 
    ps = [Ω , Γ , η , s]
    pn = [Ω_, Γ_, η_, 1/√2]
    ## initial state
    x0 = 5.0; p0 = 0.0
    u0 = [x0, p0, 5+x0*x0, im/2+x0*p0, 5+p0*p0] # σx = σpp = 5

    Random.seed!(1) 
    @named sys_fw = SDESystem(eqs)
    dict_fw = merge(Dict(unknowns(sys_fw) .=> u0), Dict(ps.=>pn))
    noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0; save_everystep=true)
    prob_fw = SDEProblem(sys_fw, dict_fw, (0, Tend); noise=noise)

    sol_fw = solve(prob_fw, EM(); dt=dt, saveat=T_saveat)

    t_W = noise.t # noise time steps
    W_W = noise.W # noise values
    N = length(t_W)-1
    dY_W = zeros(N)
    sol_x = real(sol_fw[x])
    for k in 1:N
        dW = W_W[k+1] - W_W[k]
        cd_p_c = sol_x[k]*√(2)*√(Γ_)

        dY_W[k] = dW + √η_*cd_p_c*dt
    end
    dYdt_data = [dY_W; dY_W[end]] / dt
    # dYdt(t) = dYdt_data[Int(floor(t/dt)) + 1] # TODO: why does this give a small difference?? Maybe the times do no perfectly match in the solver?
    dYdt(t) = dYdt_data[Int(round(t/dt)) + 1] 

    eqs_det = meanfield(ops, H, J; rates=R, order=2)
    function f_measure(lhs,rhs)
        term_ = √(η*Γ)*( lhs*a + a'*lhs - average(a+a')*lhs )*(Ydot(t) - √(η*Γ)*average(a+a'))
        term = cumulant_expansion(average(term_), 2)
        return rhs + term
    end
    eqs_kal = modify_equations(eqs_det, f_measure)
    eqs_kal_test = deepcopy(eqs_det)
    modify_equations!(eqs_kal_test, f_measure)
    @test isequal(eqs_kal_test[1], eqs_kal[1])
    @test isequal(eqs_kal_test[4], eqs_kal[4])

    Ydot(t) = dYdt(t)
    @named sys_kal = System(eqs_kal)
    prob_kal = ODEProblem(sys_kal, dict_fw, (0, Tend))

    sol_kal = solve(prob_kal, Euler(); dt=dt, saveat=T_saveat)
    # it = 100#8 # step from 7 to 8 changes drastically if floor() is used!! # TODO
    # abs.(sol_fw.u[it] .- sol_kal.u[it])
    @test sum(abs.(sol_fw.u[end] .- sol_kal.u[end])) < 1e-6

    # back propagation
    eqs_back_kal = meanfield(ops, -H, J; rates=R, order=2)
    eqs_back_kal_c = complete(eqs_back_kal)
    ## adjust recycling term for the back propagtion
    function f_back_lind(lhs, rhs)
        term_1 = Γ*a*lhs*a' - Γ*a'lhs*a # adapt recycling term
        term_2 = -Γ*(average(a*a') - average(a'a))*lhs
        term = cumulant_expansion(average(term_1 + term_2), 2)
        return rhs + term
    end
    eqs_back_kal_c1 = modify_equations(eqs_back_kal_c, f_back_lind)
    ## include measurement backaction 
    function f_back_kal(lhs,rhs)
        term_ = √(η*Γ)*( lhs*a' + a*lhs - average(a+a')*lhs )*(Ydot_rev(t) - √(η*Γ)*average(a+a'))
        term = cumulant_expansion(average(term_), 2)
        return rhs + term
    end
    eqs_back_kal_c2 = modify_equations(eqs_back_kal_c1, f_back_kal)
    # deterministic smoothing equations
    Ydot_rev(t) = dYdt(Tend-t)
    @named sys_back_kal = System(eqs_back_kal_c2)
    u0_bw = [0.0+0im,0,100,0,100]
    dict_bw = merge(Dict(unknowns(sys_back_kal) .=> u0_bw), Dict(ps.=>pn))
    prob_bw_kal = ODEProblem(sys_back_kal, dict_bw, (0, Tend))
    sol_bw = solve(prob_bw_kal, Euler(); dt=dt, saveat=T_saveat)

    # TODO: modify_equations!(); stochastic version - compare and test functions; 
end
