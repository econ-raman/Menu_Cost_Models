##############################################################################
##
## Type
##
##############################################################################



@with_kw struct model_params
    β::Float64 = 0.96^(1.0/12.0);        # discount factor, since will estimate at monthly level
    μ::Float64  =.002;                   # money growth rate
    σ_ϵ::Float64    =.0037;              # standard deviation of money growth rate, to match std of % of nominal GDP growth
    θ::Float64  =7;                      # Elasticity of substitution (implies markup = 7/6-1)
    ω =(θ-1)/θ;                          # Disutility of leisure, If you want 1/3 labor supply, set disutility to 3 but then you'll need lots of 3's floating around the program because log prices won't have mean 0. If you also multiply menu cost by 1/3, nothing is changed
    ρ_z = .7;                            # Persistence of productivity, estimated to match micro data levelstatistics,
    fmenu::Float64 = .045;                # Menu cost
    σ_z=.04;                              # Standard deviation of idiosyncratic productivity
end

@with_kw struct grids
    ρ_z = .7;                            # Persistence of productivity, estimated to match micro data levelstatistics,
    μ::Float64  =.002;                   # money growth rate
    npgridsize::Int64 = 181              # price gridsize, nominal price grid
    zgridsize::Int64 = 23                # idiosyncratic productivity grid size
    χgridsize::Int64 = 12;              # Aggregate state grid size (krussel-smith moment)
    ϵ_s_gridsize::Int64 = 7;            # Money shock grid size, S
    trunc=2.0;                           #Just scales grid size
    zmin::Float64 = -trunc*((1.55*.04)^2/(1-ρ_z^2))^.5;  # z is equivalent to log(z) in slides
    zmax::Float64 = trunc*((1.55*.04)^2/(1-ρ_z^2))^.5;
    npmin::Float64 = -μ*(npgridsize-1)/2;                # np is equivalent to log(p/S) in slides
    npmax::Float64 = μ*(npgridsize-1)/2;
    ϵ_s_min::Float64 = -(ϵ_s_gridsize-1)/2*μ;
    ϵ_s_max::Float64 = (ϵ_s_gridsize-1)/2*μ;
    χmin::Float64 = -.025;
    χmax::Float64 = .025;
    difftranstol::Float64 = .02;          # Krusell-Smith tolerance
    

    np::Array{Float64,1} = collect(range(npmin, npmax, length = npgridsize))
    z::Array{Float64,1} = collect(range(zmin, zmax, length = zgridsize))
    χ::Array{Float64,1} = collect(range(χmin, χmax, length = χgridsize))
    ϵ_s::Array{Float64,1} = collect(range(ϵ_s_min, ϵ_s_max, length = ϵ_s_gridsize))

end

@with_kw struct sim_param
    numfirms::Int64 = 50000
    numsim::Int64 = 500
    burnin::Int64 = 50
end
mutable struct transition_matrices
    Probz::Array{Float64,2}
    Probzcum::Array{Float64,2}
    Probϵ_s::Array{Float64,1}
    Probϵ_scum::Array{Float64,1}
end

mutable struct shocks
    moneyholder::Array{Int64,1}
    moneyshockindex::Array{Int64,1}
    zshocks::Array{Float64,2}
end

function setup_transitions(grid::grids, model::model_params)
    @unpack z, zgridsize, ϵ_s, ϵ_s_gridsize = grid
    @unpack ρ_z, σ_z, σ_ϵ = model
    w=z[2]-z[1];

    Probz = zeros(zgridsize,zgridsize)
    normcdf(x,y,z) = cdf(Normal(y,z),x)
    for j=1:zgridsize
        Probz[j,1]=normcdf((z[1]-ρ_z*z[j]+w/2)/σ_z,0,1);
        Probz[j,zgridsize]=1-normcdf((z[zgridsize]-ρ_z*z[j]-w/2)/σ_z,0,1);
        for k=2:zgridsize-1
             Probz[j,k]=normcdf((z[k]-ρ_z*z[j]+w/2)/σ_z,0,1)-normcdf((z[k]-ρ_z*z[j]-w/2)/σ_z,0,1);
        end
    end


    Probzcum = zeros(zgridsize,zgridsize)
    for j=1:zgridsize
        for k=1:zgridsize
            Probzcum[j,k]=sum(Probz[j,1:k]);
        end
    end

    #money growth is iid so no rho component
    Probϵ_s = zeros(ϵ_s_gridsize)
    w=ϵ_s[2]-ϵ_s[1];
    Probϵ_s[1]=normcdf((ϵ_s[1]+w/2)/σ_ϵ,0,1);
    Probϵ_s[ϵ_s_gridsize] =1-normcdf((ϵ_s[ϵ_s_gridsize]-w/2)/σ_ϵ,0,1);
    for j=2:ϵ_s_gridsize-1
        Probϵ_s[j]=normcdf((ϵ_s[j]+w/2)/σ_ϵ,0,1)-normcdf((ϵ_s[j]-w/2)/σ_ϵ,0,1);
    end

    Probϵ_scum = zeros(ϵ_s_gridsize)
    for j=1:ϵ_s_gridsize
        Probϵ_scum[j]=sum(Probϵ_s[1:j]);
    end

    return transition_matrices(Probz, Probzcum, Probϵ_s, Probϵ_scum)
end



function draw_shocks(grid::grids, sim::sim_param, tm::transition_matrices)
    @unpack numsim, burnin, numfirms = sim
    @unpack Probϵ_scum = tm
    @unpack ϵ_s_gridsize = grid
    moneyholder = zeros(numsim+burnin)
    moneyshockindex = zeros(numsim+burnin)
    zshocks = zeros(numsim+burnin, numfirms)
    for k = 1:numsim+burnin
            moneyrand=rand(1)[1]; # Take a random draw for the aggregate shock
            if moneyrand<Probϵ_scum[1]
                moneyholder[k]=1+(1-(ϵ_s_gridsize+1)/2);
                moneyshockindex[k]=1;
            else
                for j=1:ϵ_s_gridsize-1
                    if moneyrand>Probϵ_scum[j] && moneyrand<=Probϵ_scum[j+1]
                        moneyholder[k]=1+(j+1-(ϵ_s_gridsize+1)/2);
                        moneyshockindex[k]=j+1;  # find the index in the money grid
                    end
                end
            end

         # Draw new productivity:
         for i = 1:numfirms
         zshocks[k,i]=rand(1)[1];
         end
         
    end
        return shocks(moneyholder, moneyshockindex, zshocks)
end

