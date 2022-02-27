function update_V(grid::grids, model::model_params, tm::transition_matrices, χ̂, Vnoadjustold::Array{Float64,4}, Vadjustold::Array{Float64,3}, a0, a1)

    @unpack npgridsize, zgridsize, χgridsize, np, z, χ, ϵ_s_gridsize = grid
    @unpack ω, θ, β, fmenu = model
    @unpack Probϵ_s, Probz = tm
    fmenu = [1000, 0]
    α = [0.9, 0.1]

    Vnoadjustnew=zeros(npgridsize,zgridsize,χgridsize,2); # fNot adjusting does not depends on menu cost, so this remains the same
    Vadjustnew=zeros(zgridsize,χgridsize, 2); # f can take 2 values
    policyadjust=zeros(zgridsize,χgridsize, 2); # policy to adjust will depend on the value of f

        # Bellman of not adjustin Vn
        for i=1:npgridsize, j=1:zgridsize, k=1:χgridsize, l = 1:2
            Vnoadjustnew[i,j,k,l]=((exp(np[i]))-ω/exp(z[j]))*((exp( (a0 +a1*χ[k])*(θ-2) ))*exp(-θ*np[i]));  #Flow profit is independent of f, so this remains the same 
            #Calculate expectations/continuation values: 
            # If you don't adjust you don't care about the menu cost so this remains the same
            for n=1:zgridsize, r=1:ϵ_s_gridsize, b=1:2
                    if 1<=i-(1+(r-(ϵ_s_gridsize+1)/2)) && npgridsize>=i-(1+(r-(ϵ_s_gridsize+1)/2))  # If the money shock would push price off grid need to make sure it stays on grid
                        Vnoadjustnew[i,j,k,l]=Vnoadjustnew[i,j,k,l]+β*Probϵ_s[r]*exp(-(a0+a1*χ[k]))/exp(-(a0+a1*χ[χ̂[k,r]]))*Probz[j,n]*α[b]*max(Vnoadjustold[Int(i-1-(r-(ϵ_s_gridsize+1)/2)),n,χ̂[k,r],b],Vadjustold[n,χ̂[k,r],b]);
                    elseif 1>i-(1+(r-(ϵ_s_gridsize+1)/2))
                        Vnoadjustnew[i,j,k,l]=Vnoadjustnew[i,j,k,l]+β*Probϵ_s[r]*exp(-(a0+a1*χ[k]))/exp(-(a0+a1*χ[χ̂[k,r]]))*Probz[j,n]*α[b]*max(Vnoadjustold[1,n,χ̂[k,r],b],Vadjustold[n,χ̂[k,r],b]);
                    else
                        Vnoadjustnew[i,j,k,l]=Vnoadjustnew[i,j,k,l]+β*Probϵ_s[r]*exp(-(a0+a1*χ[k]))/exp(-(a0+a1*χ[χ̂[k,r]]))*Probz[j,n]*α[b]*max(Vnoadjustold[npgridsize,n,χ̂[k,r],b],Vadjustold[n,χ̂[k,r],b]);
                    end
            end
        end



    #The value of adjusting will be the value of not adjusting at the
    #computed at the highest value price, minus the cost of adjusting
    #to that price
    for j=1:zgridsize, k=1:χgridsize, l = 1:2
            (Vadjustnew[j,k,l], policyadjust[j,k,l])=findmax(Vnoadjustnew[:,j,k,l]);
            Vadjustnew[j,k,l]=Vadjustnew[j,k,l]-fmenu[l]*ω*exp( -(a0+a1*χ[k]) );
    end
    return Vadjustnew, Vnoadjustnew, policyadjust
end
