# Author: Raman Singh Chhina
# University of Chicago
# Purpose: To simulate the Golosov Lucas model
# February 20222

using Parameters, Distributions, Plots

include("functions.jl")
include("./VFI/update_V.jl")
include("./VFI/VFI.jl")
include("./Simulate/sim.jl")
include("./Regress/reg.jl")



# Initial guess for krusell smith:
a0=-3.3915e-4;
a1=0.5931;
r2a = 0;


using Random
Random.seed!(3)


model = model_params();
grid = grids();
sim = sim_param();
tm = setup_transitions(grid, model);
shock = draw_shocks(grid, sim, tm);

# objects that i want to see after the loop finishes
@unpack npgridsize, zgridsize, χgridsize, ϵ_s_gridsize = grid;
policy=zeros(npgridsize,zgridsize,χgridsize);
policyadjust=zeros(zgridsize,χgridsize);

χ̂ = zeros(χgridsize,ϵ_s_gridsize); # Aggregate state

# Draw and keep the shocks


# Outer loop of krussell smith: this is where we check the coefficient of
# a's. Look for a fixed point for coefficients in a linear forecasting rule

difftrans = 10
iteration = 0
#while difftrans>grid.difftranstol
while difftrans > 0.01
    iteration = iteration + 1
    println("KS iteration no", iteration)

    #Need χ̂ to stay on grid, so calculate nearest point for χ̂

    @unpack χ, zgridsize, npgridsize, χgridsize, ϵ_s_gridsize, ϵ_s = grid
    @unpack μ = model
    χ̂ = zeros(χgridsize,ϵ_s_gridsize)
    χ̂ = [a0 .+ a1 .* χ[i] - μ - ϵ_s[j] for i in 1:χgridsize, j in 1:ϵ_s_gridsize] # what would be χ̂ next period 
    χ̂ = [findmin(abs.(χ̂[i,j] .- χ))[2] for i in 1:χgridsize, j in 1:ϵ_s_gridsize] # what χ on the grid does it correspnds to

    # Solve Value Function
    Vnoadjustnew, Vadjustnew, policy = solve_policy(grid, model, tm, χ̂, a0, a1);
    
    println("Now simulating the model")
    
    # Simulate
    Csim, χsim, frequp_t, freqdown_t, sizeup_t, sizedown_t = simulate_panel(grid, sim, model, tm, χ̂, policy, shock);

    # Run Regression
    a0new, a1new, r2a = regres(sim, Csim, χsim)

    println("The r squared of the regression is ", r2a)
    # Calculate change in transition rule:  could also use absolute change instead of percentage change.
    difftrans=maximum([abs((a0-a0new)/(.5*a0+.5*a0new)) abs((a1-a1new)/(.5*a1+.5*a1new))])
    println("difftrans is ", difftrans)
    a0 = 0.25*a0 + 0.75*a0new 
    a1 = 0.25*a1 + 0.75*a1new 
end








