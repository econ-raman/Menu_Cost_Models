include("./solve_model.jl")
include("irf_sim.jl")

Csim, χsim, frequp_t, freqdown_t, sizeup_t, sizedown_t = simulate_panel(grid, sim, model, tm, χ̂, policy, shock)
Csim_ws, χsim_ws, frequp_t_ws, freqdown_t_ws, sizeup_t_ws, sizedown_t_ws = simulate_panel_extra_shock(grid, sim, model, tm, χ̂, policy, shock)


plot(Csim_ws[100:150] .- Csim[100:150])
plot(Csim_ws .- Csim)
plot(frequp_t_ws .- frequp_t)
