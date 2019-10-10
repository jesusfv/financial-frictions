% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This is the main code of the paper, which runs the loop to find the  solution of the problem

% INITIAL GUESS - by bringing it here, we make the HJB start from the previous V1 optimum every time it starts with a new K_lm
V1         = (w.*z + r.*a).^(1-gamma)/(1-gamma)/rho;  % utility from always consuming today's income
V1_stacked = V1(:);


% Firs we obtain random draws for KFE simulations (they are always the same)

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


% Creation of a grid for the PLM
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
    
    b3_HJB;                   % STEP 1: Solve the individual HJB given the PLM
    disp('Simulation started. Total number of years:')
    disp(nval_Zsim*dt)
    b5_KFE;                   % STEP 2: Simulate the path of the distribution and aggregate variables for the shocks using the KF equation
    disp('NN training started. Total number of iters:')
    disp(NN_iters*NN_starts)
    b7_PLM;                   % STEP 3: Update the PLM based on simulated data
    
    NN_starts=1;              % Only the first iteration has multiple NN starts, from the second onwards we set this to 1, as the NN uses as a first start the NN of the previous iteration

    t=0;
    % Evaluation of the error over the number of visits
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
    
    % Extra: Just to check how the NN works compared with a regression, let's run a linear regression and calculate its RMSE and R^2
    
    X_LR = [X0 X1 X2];
    beta = (X_LR'*X_LR)^-1*X_LR'*Y;
    Y_LR_fit = X_LR*beta;
    Y_LR_RMSE = (mean((Y-Y_LR_fit).^2))^0.5;
    Y_LR_R_squared = 1- ( sum((Y-Y_LR_fit).^2) ) / ( sum((Y-mean(Y)).^2) );

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

