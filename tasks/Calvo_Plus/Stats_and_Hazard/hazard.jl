



# Put them in bins
function bin_pricegaps!(stats_df::DataFrame)
    min_pg = minimum(stats_df.price_gap) 
    max_pg = maximum(stats_df.price_gap)

    pricegapbins = collect(range(min_pg, max_pg, length = 100))

    bin = zeros(size(stats_df,1))
    cpg = zeros(size(stats_df,1))

    for i in 1:size(stats_df,1)
        bin[i] = findmin(abs.(stats_df.price_gap[i] .- pricegapbins))[2]
        cpg[i] = pricegapbins[Int.(bin[i])]
    end

    stats_df.pg_bin = bin
    stats_df.pgbin_val = cpg

    return stats_df
end

function hazard(grid::grids, sim::sim_param, model::model_params, tm::transition_matrices, χ̂::Array{Int64,2}, policy::Array{Float64,4}, policyadjust::Array{Float64,3}, shock::shocks)

    @unpack numfirms, numsim, burnin = sim
    @unpack npgridsize, zgridsize, χgridsize, np, ϵ_s_gridsize, μ = grid
    @unpack Probϵ_scum, Probzcum = tm
    @unpack θ = model
    @unpack moneyholder, moneyshockindex, zshocks = shock
    Csim=zeros(numsim,1);  # Simulated demand

    #Initialize state variables
    currentfirmstate=zeros(numfirms,4);
    newfirmstate=zeros(numfirms,4);
 
        # Initial value for all the firm states is just the average value
    currentfirmstate[:,1] .= round(npgridsize/2);
    currentfirmstate[:,2] .= round(zgridsize/2);
    currentfirmstate[:,3] .= round(χgridsize/2);
    currentfirmstate[:,4] .= [ones(Int(numfirms/2)); 2 .* ones(Int(numfirms/2))]

    # Just some holders for calculating different statistics that could
    # sometimes be of interest:
 
    

    # Initialize some variables for statistics

    #change_size = zeros(numsim, numfirms,3)
    stats_df = DataFrame(time = Int[], firm = Int[], price_gap = Float64[], price_change = Float64[], price_change_ind = Int[])

    for k=1:numsim+burnin # for each period
        if k > burnin && mod(k,100) ==0
            print("Iter. Time: ", k, "  ")
        end
       
    
        pricetoday=zeros(numfirms,1);
        pricetodaywitoutmenu = zeros(numfirms)
        for i=1:numfirms
            # Price today will be the price firms choose today (normalized by M)
            pricetoday[i,1]=policy[Int(currentfirmstate[i,1]),Int(currentfirmstate[i,2]),Int(currentfirmstate[i,3]), Int(currentfirmstate[i,4])]; # Calculate each firms price today given the optimal policy
            pricetodaywitoutmenu[i] = policyadjust[Int(currentfirmstate[i,2]),Int(currentfirmstate[i,3]), Int(currentfirmstate[i,4])]
            # Calculate Hazard
            if k>burnin
                price_gap = np[Int(pricetodaywitoutmenu[i,1])]-np[Int(currentfirmstate[i,1])] # Price gap
                price_change = np[Int(pricetoday[i,1])]-np[Int(currentfirmstate[i,1])] # By how much do they actuall change price
                price_change_ind = ifelse( pricetoday[i,1] == currentfirmstate[i,1], 0,1)
                push!(stats_df, [k-burnin, i, price_gap, price_change, price_change_ind ])
            end
           
            # Price tomorrow will be price today adjusted by money shock
            newfirmstate[i,1]=pricetoday[i,1]-moneyholder[k];  # log(p/S) - (μ + ϵ_m)
            if newfirmstate[i,1]<1  # make sure it stays in grid
                newfirmstate[i,1]=1;
            elseif newfirmstate[i,1]>npgridsize
                newfirmstate[i,1]=npgridsize;
            end

            if fshocks[k,i] == 1
                newfirmstate[i,4] = 1
            else
                newfirmstate[i,4] = 2
            end
            
            # Draw new productivity:
            
            if zshocks[k,i]<Probzcum[Int(currentfirmstate[i,2]),1]
                newfirmstate[i,2]=1;
            else
                for j=1:zgridsize-1
                    if zshocks[k,i]>Probzcum[Int(currentfirmstate[i,2]),j] && zshocks[k,i]<=Probzcum[Int(currentfirmstate[i,2]),j+1]
                        newfirmstate[i,2]=j+1;
                    end
                end
            end

                
        end
        
        # We have now caluclated everything for all the firms in this period
        
        
        # Update krusell-smith state:
        newfirmstate[:,3] .= χ̂[Int(currentfirmstate[1,3]),moneyshockindex[k]]; # update the firms state # I'm not passing Xi prime to this function!!!
        if k>burnin
            Csim[k-burnin]=1/(sum(exp.(np[Int.(pricetoday[:,1])]).^(1-θ))/numfirms)^(1/(1-θ));
        end
        currentfirmstate=newfirmstate;
        
    end 
    return stats_df

end