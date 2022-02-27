# Function to solve the final regression for Krussel smith


function regres(sim::sim_param, Csim::Array{Float64}, χsim::Array{Float64})
    # Setup matrix for krusell-smith regression
    @unpack numsim, numfirms = sim
    logPoverMregression = zeros(numsim,3)
    for k=1:numsim
        logPoverMregression[k,1]=log(1/Csim[k]);  #y variable
        logPoverMregression[k,2]=1;  #constant
        logPoverMregression[k,3]= χsim[k];  #x variable

    end
    logPoverMregression = logPoverMregression[1:end,:]
    # (X'X)^(-1)*X'Y:
    coefficients=inv(logPoverMregression[:,2:3]'*logPoverMregression[:,2:3])*(logPoverMregression[:,2:3]'*logPoverMregression[:,1]);
    a0new=coefficients[1];
    a1new=coefficients[2];
    yhat=logPoverMregression[:,2:3]*coefficients;
    r2a=1-sum((logPoverMregression[:,1]-yhat).^2)/sum((logPoverMregression[:,1] .- sum(logPoverMregression[:,1]/length(logPoverMregression[:,1]))).^2)
    return a0new, a1new, r2a
end