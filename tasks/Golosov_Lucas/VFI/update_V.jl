
function update_V(grid::grids, model::model_params, tm::transition_matrices, χ̂, Vnoadjustold::Array{Float64,3}, Vadjustold::Array{Float64,2}, a0, a1)

    @unpack npgridsize, zgridsize, χgridsize, np, z, χ, ϵ_s_gridsize = grid
    @unpack ω, θ, β, fmenu = model
    @unpack Probϵ_s, Probz = tm

    Vnoadjustnew=zeros(npgridsize,zgridsize,χgridsize);
    Vadjustnew=zeros(zgridsize,χgridsize);
    policyadjust=zeros(zgridsize,χgridsize);

        # Bellman of not adjustin Vn
        for i=1:npgridsize, j=1:zgridsize, k=1:χgridsize
            Vnoadjustnew[i,j,k]=((exp(np[i]))-ω/exp(z[j]))*((exp( (a0 +a1*χ[k])*(θ-2) ))*exp(-θ*np[i]));  #Flow profit
            #Calculate expectations/continuation values:
            for n=1:zgridsize, r=1:ϵ_s_gridsize
                    if 1<=i-(1+(r-(ϵ_s_gridsize+1)/2)) && npgridsize>=i-(1+(r-(ϵ_s_gridsize+1)/2))  # If the money shock would push price off grid need to make sure it stays on grid
                        Vnoadjustnew[i,j,k]=Vnoadjustnew[i,j,k]+β*Probϵ_s[r]*exp(-(a0+a1*χ[k]))/exp(-(a0+a1*χ[χ̂[k,r]]))*Probz[j,n]*max(Vnoadjustold[Int(i-1-(r-(ϵ_s_gridsize+1)/2)),n,χ̂[k,r]],Vadjustold[n,χ̂[k,r]]);
                    elseif 1>i-(1+(r-(ϵ_s_gridsize+1)/2))
                        Vnoadjustnew[i,j,k]=Vnoadjustnew[i,j,k]+β*Probϵ_s[r]*exp(-(a0+a1*χ[k]))/exp(-(a0+a1*χ[χ̂[k,r]]))*Probz[j,n]*max(Vnoadjustold[1,n,χ̂[k,r]],Vadjustold[n,χ̂[k,r]]);
                    else
                        Vnoadjustnew[i,j,k]=Vnoadjustnew[i,j,k]+β*Probϵ_s[r]*exp(-(a0+a1*χ[k]))/exp(-(a0+a1*χ[χ̂[k,r]]))*Probz[j,n]*max(Vnoadjustold[npgridsize,n,χ̂[k,r]],Vadjustold[n,χ̂[k,r]]);
                    end
            end
        end



    #The value of adjusting will be the value of not adjusting at the
    #computed at the highest value price, minus the cost of adjusting
    #to that price
    for j=1:zgridsize     
        for k=1:χgridsize
            (Vadjustnew[j,k], policyadjust[j,k])=findmax(Vnoadjustnew[:,j,k]);
            Vadjustnew[j,k]=Vadjustnew[j,k]-fmenu*ω*exp( -(a0+a1*χ[k]) );
        end
    end
    return Vadjustnew, Vnoadjustnew, policyadjust
end