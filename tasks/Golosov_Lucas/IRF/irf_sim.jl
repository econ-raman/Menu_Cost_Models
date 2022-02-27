
function simulate_panel_extra_shock(grid::grids, sim::sim_param, model::model_params, tm::transition_matrices, χ̂::Array{Int64,2}, policy::Array{Float64,3}, shock::shocks)

    @unpack numfirms, numsim, burnin = sim
    @unpack npgridsize, zgridsize, χgridsize, np, ϵ_s_gridsize, μ = grid
    @unpack Probϵ_scum, Probzcum = tm
    @unpack θ = model
    @unpack moneyholder, moneyshockindex, zshocks = shock
    #Initialize state variables
    currentfirmstate=zeros(numfirms,3);
    newfirmstate=zeros(numfirms,3);

        # Initial value for all the firm states is just the average value
    currentfirmstate[:,1] .= round(npgridsize/2);
    currentfirmstate[:,2] .= round(zgridsize/2);
    currentfirmstate[:,3] .= round(χgridsize/2);

    numup=0;
    numdown=0;
    sizeup=0;
    sizedown=0;

    Csim=zeros(numsim,1);  # Simulated demand
    χsim=zeros(numsim,1);  # Simulated Krusell-Smith moment
    eps=zeros(numsim,1); # Money shock

    # Just some holders for calculating different statistics that could
    # sometimes be of interest:
    frequp_t=zeros(numsim,1);
    freqdown_t=zeros(numsim,1);
    sizeup_t=zeros(numsim,1);
    sizedown_t=zeros(numsim,1);
    size_t=zeros(numsim,1);

    inflation=zeros(numsim-1,1);  # keeps track of implied inflation

    # Initialize some variables for statistics
    xsd=0;
    xsdholder=zeros(numsim,1);
    meanpricechange=zeros(numsim,1);

    for k=1:numsim+burnin # for each period

        if k > burnin && mod(k,100) ==0
            print("Iter. Time: ", k, "  ")
        end 

        numup2=0;
        numdown2=0;
        sizeup2=0;
        sizedown2=0;
        size2=0;

        if k == burnin + 100
            if moneyshockindex[k] == ϵ_s_gridsize
                error("the money shock already had the highest value ")
            end
            moneyshockindex[k] = moneyshockindex[k] + 1;
            moneyholder[k] = moneyholder[k] + 1
        end
        
      
        pricetoday=zeros(numfirms,1);
        for i=1:numfirms
            # Price today will be the price firms choose today (normalized by M)
            pricetoday[i,1]=policy[Int(currentfirmstate[i,1]),Int(currentfirmstate[i,2]),Int(currentfirmstate[i,3])]; # Calculate each firms price today given the optimal policy
            if k>burnin  # Calculate statistics:
                if pricetoday[i,1]>currentfirmstate[i,1] # The firms which decrease the price
                    numup=numup+1;
                    numup2=numup2+1;
                    sizeup=sizeup+np[Int(pricetoday[i,1])]-np[Int(currentfirmstate[i,1])];
                    sizeup2=sizeup2+np[Int(pricetoday[i,1])]-np[Int(currentfirmstate[i,1])];
                    size2=size2+np[Int(pricetoday[i,1])]-np[Int(currentfirmstate[i,1])];
                elseif pricetoday[i,1]<currentfirmstate[i,1] # The firms which increase the price
                    numdown=numdown+1;
                    numdown2=numdown2+1;
                    sizedown=sizedown+np[Int(currentfirmstate[i,1])]-np[Int(pricetoday[i,1])];
                    sizedown2=sizedown2+np[Int(currentfirmstate[i,1])]-np[Int(pricetoday[i,1])];
                    size2=size2+np[Int(pricetoday[i,1])]-np[Int(currentfirmstate[i,1])];
                end
            end

            
            
            # Price tomorrow will be price today adjusted by money shock
            newfirmstate[i,1]=pricetoday[i,1]-moneyholder[k];  # log(p/S) - (μ + ϵ_m)
            if newfirmstate[i,1]<1  # make sure it stays in grid
                newfirmstate[i,1]=1;
            elseif newfirmstate[i,1]>npgridsize
                newfirmstate[i,1]=npgridsize;
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
        # calculate stats:
        if k>burnin
            xsd=xsd+sum(((np[Int.(pricetoday[:,1])] .- np[Int.(currentfirmstate[:,1])]) .- mean(np[Int.(pricetoday[:,1])] .- np[Int.(currentfirmstate[:,1])])).^2/numfirms)^.5;
            xsdholder[k-burnin]=sum(((np[Int.(pricetoday[:,1])] .- np[Int.(currentfirmstate[:,1])]) .- mean(np[Int.(pricetoday[:,1])] .- np[Int.(currentfirmstate[:,1])])).^2/numfirms)^.5;
            meanpricechange[k-burnin]=mean(np[Int.(pricetoday[:,1])] .- np[Int.(currentfirmstate[:,1])]);
        end
        
        # Update krusell-smith state:
        newfirmstate[:,3] .= χ̂[Int(currentfirmstate[1,3]),moneyshockindex[k]]; # update the firms state # I'm not passing Xi prime to this function!!!
        
        # calculate simulated variables
        if k>burnin
            Csim[k-burnin]=1/(sum(exp.(np[Int.(pricetoday[:,1])]).^(1-θ))/numfirms)^(1/(1-θ));
            χsim[k-burnin,1]=log((sum((exp.(np[Int.(currentfirmstate[:,1])])).^(1-θ))/numfirms)^(1/(1-θ)));
            eps[k-burnin,1]=moneyholder[k]* μ;
            frequp_t[k-burnin,1]=numup2/numfirms;
            freqdown_t[k-burnin,1]=numdown2/numfirms;
            sizeup_t[k-burnin,1]=sizeup2/numup2;
            sizedown_t[k-burnin,1]=sizedown2/numdown2;
            size_t[k-burnin,1]=size2/(numup2+numdown2);
        end
        if k>burnin+1
            changemv=eps[k-burnin-1,1];
            changec=log(Csim[k-burnin]/Csim[k-burnin-1]);
            inflation[k-burnin-1]=changemv-changec;
        end

        
        
        currentfirmstate=copy(newfirmstate);
        
    end
    return Csim, χsim, frequp_t, freqdown_t, sizeup_t, sizedown_t

end

