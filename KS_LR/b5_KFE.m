% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This code simulates the complete path of the distribution given the
% shock


Ksim    = zeros(nval_Zsim,1);
KposU   = zeros(nval_Zsim,1);
KposD   = zeros(nval_Zsim,1);

KKposU  = zeros(nval_Zsim,1);
KKposD  = zeros(nval_Zsim,1);

rsim    = zeros(nval_Zsim,1);

g_big   = zeros(nval_a,nval_z,nval_Zsim);

% MAIN LOOP

disp('KFE computations will NOT use the parallel toolbox');

for t=1:nval_Zsim

    if t==1
        g1=g_ss;
    else

        g0=g1;
        
        % Interpolation step to compute next period's distribution
        % The following block computes the transition g(t+1) = (I - A*dt)^(-1)*g(t) for the four closest grid points

        myA        = A1{KposD(t-1,1),ZposD(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1DD       = reshape(g1_stacked,size(g0));

        myA        = A1{KposD(t-1,1),ZposU(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1DU       = reshape(g1_stacked,size(g0));

        myA        = A1{KposU(t-1,1),ZposD(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1UD       = reshape(g1_stacked,size(g0));

        myA        = A1{KposU(t-1,1),ZposU(t-1,1)};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1UU       = reshape(g1_stacked,size(g0));

        % The we compute distribution next period, averaging those four results
        g1 = wK*wZ(t-1)*g1DD + wK*(1-wZ(t-1))*g1DU + (1-wK)*wZ(t-1)*g1UD + (1-wK)*(1-wZ(t-1))*g1UU;

    end
    
    g_big(:,:,t)=g1;
    
    Ksim(t)    =sum(sum(g1.*a2*da));
    
    Ksim(t)  =   max([Ksim(t) Kmin+0.000001]);
    Ksim(t)  =   min([Ksim(t) Kmax-0.000001]);
    KposD(t) = floor((Ksim(t)-Kmin)/dK)+1;
    KposU(t) =  ceil((Ksim(t)-Kmin)/dK)+1;
    wK=(K_grid(KposU(t))-Ksim(t))/dK;         % weight of KposD

    KKposD(t) = floor((Ksim(t)-Kmin)/dKK)+1;
    KKposU(t) =  ceil((Ksim(t)-Kmin)/dKK)+1;
    
end



