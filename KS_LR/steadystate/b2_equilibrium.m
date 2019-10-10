function K = b2_equilibrium(parameters,K)
%Computes the excess demand
maxitK  = parameters.maxitK;
critK   = parameters.critK;
weK     = parameters.weK; 

for n=1:maxitK
    results      = b3_HJB_stationary(parameters,K);
    distribution = b4_KFE_stationary(parameters,results);
    S = sum(sum(results.aa .*distribution *results.da));
    if abs(S-K)<critK
        break
    else
        K = weK * S +(1-weK)* K;
        disp(K)
    end
end

if n==maxitK
    disp('ERROR: max iterations reached')
end
