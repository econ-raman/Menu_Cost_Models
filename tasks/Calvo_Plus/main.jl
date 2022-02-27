using Parameters, Distributions, Plots, DataFrames

include("solve_model.jl")

using Random
Random.seed!(3)


model = model_params();
grid = grids();
sim = sim_param();
tm = setup_transitions(grid, model);
shock = draw_shocks(grid, sim, tm);

# Draw and keep the menu cost shocks

fshocks = zeros(sim.numsim+sim.burnin, sim.numfirms)
for k in 1:sim.numsim+sim.burnin-1, i = 1:sim.numfirms
    j=rand(1)[1];
    if j < 0.9
        fshocks[k,i] = 1
    else
        fshocks[k,i] = 2
    end
end

# Solve the model
policy, policyadjust, a0, a1, χ̂ = solve_model(a0_initial, a1_initial, model, grid, shock, tm, sim)

# ------------------------------------------ Compute IRF --------------------------------------------------------- #

include("./IRF/irf_sim.jl")

Csim, χsim, frequp_t, freqdown_t, sizeup_t, sizedown_t = simulate_panel(grid, sim, model, tm, χ̂, policy, shock)
Csim_ws, χsim_ws, frequp_t_ws, freqdown_t_ws, sizeup_t_ws, sizedown_t_ws = simulate_panel_extra_shock(grid, sim, model, tm, χ̂, policy, shock)

# Plot the IRF
plot(Csim_ws[100:150] .- Csim[100:150], linewidth = 3, label = "", xlabel = "Time", ylabel = "Impulse Response")
savefig("./output/C_irf_50_periods.png")

# Plot the IRF
plot(Csim_ws .- Csim, linewidth = 3, label = "", xlabel = "Time", ylabel = "Impulse Response")
savefig("./output/C_irf_500_periods.png")

# Average IRF

avg_IRF = mean( log.(Csim_ws[100]/Csim[100])) / (grid.ϵ_s[2] - grid.ϵ_s[1])

# Variance

var_C = var(Csim)

# Double the menu cost 
model2 = model_params(fmenu = 2 * model.fmenu)
policy, policyadjust, a0, a1, χ̂ = solve_model(a0_initial, a1_initial, model2, grid, shock, tm, sim) # Solve the model again
Csim, χsim, frequp_t, freqdown_t, sizeup_t, sizedown_t = simulate_panel(grid, sim, model, tm, χ̂, policy, shock)

var_C2 = var(Csim)

# ---------------------------------- Hazard ------------------------------------------ #

include("./Stats_and_Hazard/hazard.jl")

stats_df = hazard(grid, sim, model, tm, χ̂, policy, policyadjust, shock)
bin_pricegaps!(stats_df)

hazard_df = combine(groupby(stats_df, :pg_bin), :price_change_ind => mean, :pgbin_val => mean)
scatter(hazard_df.pgbin_val_mean, hazard_df.price_change_ind_mean, label = "", markershape = :hexagon) # Plot hazard
xlabel!("price gap")
ylabel!("hazard")

savefig("./output/hazard_cp.png")

# ---------------------------- Distribution of observed price changes --------------------------- #

dist_prices = combine(groupby(stats_df, :price_change), nrow)
sort!(dist_prices)

dist_prices_without_zero = dist_prices[dist_prices.price_change .!== 0.0,:]
bar(dist_prices_without_zero.price_change, dist_prices_without_zero.x1, label = "")
xlabel!("observed price changes")
ylabel!("frequency")
savefig("./output/observed_p_changes_cp.png")

# ------------------------------ Distribution of price gaps -------------------------------------- #

price_gap_dist = combine(groupby(stats_df, :pgbin_val), nrow)
bar(price_gap_dist.pgbin_val, price_gap_dist.x1 ./ (sum(price_gap_dist.x1)), label = "")
xlabel!("price gap")
ylabel!("density")
savefig("./output/price_gap_dist_cp.png")

