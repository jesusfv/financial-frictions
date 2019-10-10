% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This is the main code of the paper, which runs the loop to find the  solution of the problem

% INITIAL GUESS - by bringing it here, we make the HJB start from the previous V1 optimum every time it starts with a new K_lm
V1         = (w.*z + r.*a).^(1-gamma)/(1-gamma)/rho;  % utility from always consuming today's income
V1_stacked = V1(:);


% random draws for KFE simulations (they are always the same)

Zshocks = randn(nval_Zsim,1);

Zsim    = zeros(nval_Zsim,1);
ZposU   = zeros(nval_Zsim,1);
ZposD   = zeros(nval_Zsim,1);
wZ      = zeros(nval_Zsim,1);

Zsim(1,1)=Zmean;
for t=2:nval_Zsim
    Zsim(t,1)=theta*dt*Zmean+(1-theta*dt)*Zsim(t-1,1)+sigma*sqrt(dt)*Zshocks(t);
end

for t=1:nval_Zsim
    Zsim(t)  =   max([Zsim(t) Zmin+0.000001]);
    Zsim(t)  =   min([Zsim(t) Zmax-0.000001]);
    ZposD(t) = floor((Zsim(t)-Zmin)/dZ)+1;
    ZposU(t) =  ceil((Zsim(t)-Zmin)/dZ)+1;
    wZ(t)    = (Z_grid(ZposU(t))-Zsim(t))/dZ;    % weight of ZposD
end


% PLM iteration
PLM   = zeros(nval_a, nval_z, nval_K, nval_Z);   % if left at zero, the first guess is "no expected change in K", for all points
PLM_2 = zeros(nval_a, nval_z, nval_K, nval_Z);   % (the PLM will always be flat on its first two dimmensions; but it's easier to have it as 4D)

PLM_finegrid   = zeros(nval_KK, nval_Z);         % finer grid, used for training the NN, for determining visited range and for the convergence criteria
PLM_finegrid_2 = zeros(nval_KK, nval_Z);
PLM_visits     = zeros(nval_KK, nval_Z).*NaN;    % Number of visits in each fine-grid point in the current iteration


myfig=figure(1);
surf(Z_grid,KK_grid,PLM_finegrid)
title('Perceived Law of Motion','interpreter','Latex')
xlabel('Z','interpreter','Latex')
ylabel('K','interpreter','Latex')
axis([ Zmin Zmax Kmin Kmax -0.15 0.15])
savefig(myfig,'fig0.fig');
close all

X0 = ones(nval_Zsim-delay_Zsim,1);      % we will use this as constant in the just-for-reporting linear regression

format long;
diff_vector = [];


% MAIN LOOP: this loop iterates until maxitPLM or convergence in the PLM
for it_PLM=1:maxitPLM

    PLM_visits_=PLM_visits;   % Number of visits in each grid point in the previous iteration
    
    b3_HJB;                   % Solve the individual HJB given the PLM
    disp('Simulation started. Total number of years:')
    disp(nval_Zsim*dt)
    b5_KFE;                   % Simulate the path of the distribution and aggregate variables for the shocks
    
    % PLM: linear regression
    
    Y = (Ksim(delay_Zsim+1:nval_Zsim)-Ksim(delay_Zsim:nval_Zsim-1))/dt; % dK(t): Growth rate of K at time t
    X1=  Ksim(delay_Zsim:nval_Zsim-1);                                  % K(t-1)
    X2=  Zsim(delay_Zsim:nval_Zsim-1);                                  % Z(t-1)
    
    X1=log(X1);
    X3=X1.*X2;

    X_LR = [X0 X1 X2 X3];
    beta = (X_LR'*X_LR)^-1*X_LR'*Y;
    Y_LR_fit = X_LR*beta;
    Y_LR_RMSE = (mean((Y-Y_LR_fit).^2))^0.5;
    Y_LR_R_squared = 1- ( sum((Y-Y_LR_fit).^2) ) / ( sum((Y-mean(Y)).^2) );
    
    % copy all this into all other items, so I can reuse the code from the NN implementation without any further changes
    
    Y_fit_IP   = Y_LR_fit;
    Y_error_IP = Y-Y_fit_IP;
    Y_RMSE_IP  = (mean(Y_error_IP.^2))^0.5;
    
    Y_RMSE_NN1 = Y_LR_RMSE;
    Y_fit      = Y_LR_fit;
    Y_RMSE_NN2 = Y_LR_RMSE;
    Y_R_squared= 1- ( sum((Y-Y_fit).^2) ) / ( sum((Y-mean(Y)).^2) ); % This computes the R^2 over the whole sample (all the points in the simulation, except the discarded ones at the beginning)

    % Use LR to calculate PLM, first on fine grid then on HJB grid

    X1mm=    KK_grid_2D ;
    X2mm=    Z_grid_2D ;

    X1m = reshape(X1mm,[nval_KK*nval_Z,1]);
    X2m = reshape(X2mm,[nval_KK*nval_Z,1]);
    X0m = ones(size(X1m));
    
    X1m=log(X1m);
    X3m=X1m.*X2m;

    X_LR = [X0m X1m X2m X3m];
    Y_LR_fit = X_LR*beta;
    PLM_finegrid_2 = reshape(Y_LR_fit,size(X1mm));

    X1mm=    squeeze(K(1,1,:,:)) ;
    X2mm=    squeeze(Z(1,1,:,:)) ;

    X1m = reshape(X1mm,[nval_K*nval_Z,1]);
    X2m = reshape(X2mm,[nval_K*nval_Z,1]);
    X0m = ones(size(X1m));
    
    X1m=log(X1m);
    X3m=X1m.*X2m;

    X_LR = [X0m X1m X2m X3m];
    Y_LR_fit = X_LR*beta;
    Y_mat2s  = reshape(Y_LR_fit,size(X1mm));

    for ia=1:nval_a
        for iz=1:nval_z
            PLM_2(ia,iz,:,:)=Y_mat2s(:,:);
        end
    end

    PLM_IP = PLM_2;
    
    % continue

    t=0;
    PLM_visits=zeros(nval_KK,nval_Z);      % Reset the number of visits: Let's count how many times we visited each point in the fine grid
    for t=delay_Zsim+1:nval_Zsim           % Looping over the used simulation points, we will count "0.25 visits" for each of the four closest grid points to each simulated point
        PLM_visits(KKposD(t-1),ZposU(t-1))=PLM_visits(KKposD(t-1),ZposU(t-1))+0.25;
        PLM_visits(KKposU(t-1),ZposU(t-1))=PLM_visits(KKposU(t-1),ZposU(t-1))+0.25;
        PLM_visits(KKposD(t-1),ZposD(t-1))=PLM_visits(KKposD(t-1),ZposD(t-1))+0.25;
        PLM_visits(KKposU(t-1),ZposD(t-1))=PLM_visits(KKposU(t-1),ZposD(t-1))+0.25;
    end
    mysum_ = sum(sum(PLM_visits));       % Total number of visits per grid point: We use this one for the RMSdiff, which is weighted by number of visits to each point

    mysum  = nval_Z*nval_KK;            % mysum is the number of fine-grid points that have at least one visit; we start assuming all of them are visited
    PLM_diff=PLM_finegrid_2-PLM_finegrid;
    for iZ=1:nval_Z                                
        for iK=1:nval_KK
            if PLM_visits(iK,iZ)==0      % if this grid point has not been visited in the current set of simulations...
                if PLM_visits_(iK,iZ)==0 % ... and it also was not visited in the previous set of simulations...
                    PLM_diff(iK,iZ)=0;   % ... then this point is ignored for diff criteria...
                    mysum=mysum-1;       % ... and the number of visited points falls by one
                end
           end
        end
    end
    diff_uw = (sum(sum( PLM_diff             .^2))/mysum )^0.5;      %RMSdiff computed only over the visited points of the fine grid - unweighted
    diff    = (sum(sum((PLM_diff.*PLM_visits).^2))/mysum_)^0.5;      %RMSdiff computed only over the visited points of the fine grid - weighted by number of visits to each point in that grid
    diff_vector = [diff_vector;diff];
    
    % This is just for the plotting routines
    for iZ=1:nval_Z
        for iK=1:nval_KK
            if PLM_visits(iK,iZ)==0
                PLM_visits(iK,iZ)=NaN;   % substitute zeros for NaN so it looks better in the graphs
            end
        end
    end
    
    % report summary results from this iteration
    disp('   iter                Y_LR_RMSE           Y_RMSE_NN2           Y_RMSE_IP           diff_weighted:')
    disp([   it_PLM              Y_LR_RMSE           Y_RMSE_NN2           Y_RMSE_IP           diff          ]);

    % UPDATE PLM using a relaxation algorithm
    PLM          = (1-wePLM)*PLM          + wePLM*PLM_2;
    PLM_finegrid = (1-wePLM)*PLM_finegrid + wePLM*PLM_finegrid_2;

    myfig=figure(1);
    surf(Z_grid,KK_grid,PLM_finegrid)
    title('Perceived Law of Motion','interpreter','Latex')
    xlabel('Z','interpreter','Latex')
    ylabel('K','interpreter','Latex')
    axis([ Zmin Zmax Kmin Kmax -0.15 0.15])
    savefig(myfig,strcat('fig',num2str(it_PLM),'.fig'));
    close all
    
    for iZ=1:nval_Z
        for iK=1:nval_KK
           if and(isnan(PLM_visits(iK,iZ)),isnan(PLM_visits_(iK,iZ)))
               PLM_diff(iK,iZ)=NaN;                                   % put NaN on non-visited points, so they don't appear in the graph
           end
        end
    end
    
    myfig=figure(1);
    surf(Z_grid,KK_grid,PLM_diff);
    saveas(myfig,strcat('figdif',num2str(it_PLM),'.fig'))
    close all
    
    wePLM = wePLM*wePLM1+wePLM2;
    
    if and(diff<critPLM,it_PLM>5)  % stop if diff meets criterion, but not if we've ran surprisingly few iterations

        PLM = PLM_2;
        PLM_finegrid = PLM_finegrid_2;

        myfig=figure(1);
        surf(Z_grid,KK_grid,PLM_finegrid)
        title('Perceived Law of Motion','interpreter','Latex')
        xlabel('Z','interpreter','Latex')
        ylabel('K','interpreter','Latex')
        axis([ Zmin Zmax Kmin Kmax -0.15 0.15])
        savefig(myfig,strcat('fig',num2str(it_PLM+1),'.fig'));
        close all

        b3_HJB;
        b5_KFE;

        break

    end
end

if it_PLM==maxitPLM
    disp('ERROR: max iterations reached')
end

