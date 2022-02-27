function solve_policy(grid::grids, model::model_params, tm::transition_matrices, χ̂, a0, a1)
    @unpack npgridsize, zgridsize, χgridsize, np, z, χ, ϵ_s_gridsize = grid

    #Initialize guesses and matrices:
    

    Vadjustold=zeros(zgridsize,χgridsize,2);
    policyadjust=zeros(zgridsize,χgridsize,2);
    Vnoadjustold=zeros(npgridsize,zgridsize,χgridsize,2);
    policy=zeros(npgridsize,zgridsize,χgridsize,2);
    policyold=copy(policy);

    # Holder for V new'S

    Vnoadjustnew=zeros(npgridsize,zgridsize,χgridsize,2); # the states are nominal price level of the firm, idiosyncratic shock z and χ
        Vadjustnew=zeros(zgridsize,χgridsize,2);
    
    diffmax=10;
    
    iter = 0
    while diffmax>0  
        iter = iter + 1
            

        Vadjustnew, Vnoadjustnew, policyadjust = update_V(grid, model, tm, χ̂, Vnoadjustold, Vadjustold, a0, a1) # Given the coefficients a0 and a1 we update the value function
    
      #Policy is max of adjusting and not adjusting:
      for i=1:npgridsize, j=1:zgridsize, k=1:χgridsize, l = 1:2
        if Vnoadjustnew[i,j,k,l]>Vadjustnew[j,k,l]
            policy[i,j,k,l]=i; # if the value function for not adjusting is higher then the price should be kept the same irrespective of the productivity and χ
        else
            policy[i,j,k,l]=policyadjust[j,k,l]; # other wise the policy is to move to the optimal value
        end
    end

        
        #Converged when policy function doesn't change:
        diffmax=maximum((abs.(policy-policyold)))
        policyold= copy(policy);
        Vnoadjustold=copy(Vnoadjustnew);
        Vadjustold=copy(Vadjustnew);
    end

    println("VFI took ", iter, " iterations")
    return Vnoadjustnew, Vadjustnew, policy, policyadjust
end