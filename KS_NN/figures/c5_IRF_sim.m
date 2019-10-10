% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

irf_Zsim(1,1)=theta*dt*Zmean+(1-theta*dt)*Zmean+sigma*sqrt(dt)*e_irf(1);

for t=2:nval_IRF
    irf_Zsim(t,1)=theta*dt*Zmean+(1-theta*dt)*irf_Zsim(t-1,1)+sigma*sqrt(dt)*e_irf(t);
end

for t=1:nval_IRF
    irf_Zsim(t)  =   max([irf_Zsim(t) Zmin+0.000001]);
    irf_Zsim(t)  =   min([irf_Zsim(t) Zmax-0.000001]);
    irf_ZposD(t) = floor((irf_Zsim(t)-Zmin)/dZ)+1;
    irf_ZposU(t) =  ceil((irf_Zsim(t)-Zmin)/dZ)+1;
    irf_wZ(t)    = (Z_grid(irf_ZposU(t))-irf_Zsim(t))/dZ;    % weight of ZposD
end


for t=1:nval_IRF
    
    if t>1
        
        g0=g1;
        
        % Interpolation step to compute next period's distribution
        % The following block computes the transition g(t+1) = (I - A*dt)^(-1)*g(t) for the four closest grid points

        myA        = A1{irf_KposD(t-1,1),irf_ZposD(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1DD       = reshape(g1_stacked,size(g0));

        myA        = A1{irf_KposD(t-1,1),irf_ZposU(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1DU       = reshape(g1_stacked,size(g0));

        myA        = A1{irf_KposU(t-1,1),irf_ZposD(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1UD       = reshape(g1_stacked,size(g0));

        myA        = A1{irf_KposU(t-1,1),irf_ZposU(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1UU       = reshape(g1_stacked,size(g0));

        % The we compute distribution next period, averaging those four results
        g1 = irf_wK*irf_wZ(t-1)*g1DD + irf_wK*(1-irf_wZ(t-1))*g1DU + (1-irf_wK)*irf_wZ(t-1)*g1UD + (1-irf_wK)*(1-irf_wZ(t-1))*g1UU;

    end
    
    irf_gsim(:,:,t)=g1;
    
    irf_Ksim(t)    =sum(sum(g1.*a2*da));
    
    irf_Ksim(t)  =   max([irf_Ksim(t) Kmin+0.000001]);
    irf_Ksim(t)  =   min([irf_Ksim(t) Kmax-0.000001]);
    irf_KposD(t) = floor((irf_Ksim(t)-Kmin)/dK)+1;
    irf_KposU(t) =  ceil((irf_Ksim(t)-Kmin)/dK)+1;
    irf_wK=(K_grid(irf_KposU(t))-irf_Ksim(t))/dK;         % weight of KposD

end
