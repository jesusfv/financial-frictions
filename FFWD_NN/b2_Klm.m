% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This is the main code of the paper, which runs the loop to find the  solution of the problem

% INITIAL GUESS - by bringing it here, we make the HJB start from the previous V1 optimum every time it starts with a new PLM
V1      = ((w.*z + r.*a).^(1-gamma)-1)/(1-gamma)/rho;  % utility from always consuming today's income
% construct initial V guess in shape of V1_stacked
V1_stacked = V1(:);

% random draws for KFE simulations (they are always the same)

rng(rngseed1);
e_shocks = randn(nval_sim,1);
sdte_shocks = sigma*(dt^0.5)*e_shocks;

% PLM iteration
PLM   = zeros(nval_a, nval_z, nval_B, nval_N);   % if left at zero, the first guess is "no expected change in B", for all points
PLM_2 = zeros(nval_a, nval_z, nval_B, nval_N);   % (the PLM will always be flat on its first two dimmensions; but it's easier to have it as 4D)

PLM_finegrid   = zeros(nval_BB, nval_NN);        % finer grid, used for training the NN, for determining visited range and for the convergence criteria
PLM_finegrid_2 = zeros(nval_BB, nval_NN);
PLM_visits     = zeros(nval_BB, nval_NN).*NaN;   % Number of visits in each fine-grid point in the current iteration

% jumpstart: override first guess of the PLM (all zeros) with the result from the version of the  model with a linear regression
eval(load_PLM)

myfig=figure(1);
surf(NN_grid,BB_grid,PLM_finegrid)
title('Perceived Law of Motion','interpreter','Latex')
xlabel('N','interpreter','Latex')
ylabel('B','interpreter','Latex')
axis([ Nmin Nmax Bmin Bmax -0.15 0.15])
savefig(myfig,'fig0.fig');
close all

X0=    ones(multi_sim*used_sim,1);               % we will use this as constant in the just-for-reporting linear regression

format long;
diff_vector = [];


% MAIN LOOP: this loop iterates until maxitPLM or convergence in the PLM
for it_PLM=1:maxitPLM

    PLM_visits_=PLM_visits;   % Number of visits in each grid point in the previous iteration
    
    b3_HJB;                   % Solve the individual HJB given the PLM
    disp('Simulation started. Total number of years:')
    disp(nval_sim*dt)
    b5_KFE;                   % Simulate the path of the distribution and aggregate variables for the shocks
    disp('NN training started. Total number of iters:')
    disp(NN_iters*NN_starts)
    b7_PLM;                   % Recompute the PLM based on simulated data
    
    NN_starts=1;              % Only the first iteration has multiple NN starts, from the second onwards we set this to 1, as the NN uses as a first start the NN of the previous iteration

    t=0;
    PLM_visits=zeros(nval_BB,nval_NN);     % Reset the number of visits: Let's count how many times we visited each point in the fine grid
    for it_sim=1:multi_sim
        t=t+delay_sim;
        for it_t=1:used_sim                % Looping over the used simulation points, we will count "0.25 visits" for each of the four closest grid points to each simulated point
            t=t+1;
            PLM_visits(BBposD(t-1),NNposU(t-1))=PLM_visits(BBposD(t-1),NNposU(t-1))+0.25;
            PLM_visits(BBposU(t-1),NNposU(t-1))=PLM_visits(BBposU(t-1),NNposU(t-1))+0.25;
            PLM_visits(BBposD(t-1),NNposD(t-1))=PLM_visits(BBposD(t-1),NNposD(t-1))+0.25;
            PLM_visits(BBposU(t-1),NNposD(t-1))=PLM_visits(BBposU(t-1),NNposD(t-1))+0.25;
        end
    end
    mysum_ = sum(sum(PLM_visits));       % Total number of visits per grid point: We use this one for the RMSdiff, which is weighted by number of visits to each point

    mysum  = nval_NN*nval_BB;            % mysum is the number of fine-grid points that have at least one visit; we start assuming all of them are visited
    PLM_diff=PLM_finegrid_2-PLM_finegrid;
    for iN=1:nval_NN                                
        for iB=1:nval_BB
            if PLM_visits(iB,iN)==0      % if this grid point has not been visited in the current set of simulations...
                if PLM_visits_(iB,iN)==0 % ... and it also was not visited in the previous set of simulations...
                    PLM_diff(iB,iN)=0;   % ... then this point is ignored for diff criteria...
                    mysum=mysum-1;       % ... and the number of visited points falls by one
                end
           end
        end
    end
    diff_uw = (sum(sum( PLM_diff             .^2))/mysum )^0.5;      %RMSdiff computed only over the visited points of the fine grid - unweighted
    diff    = (sum(sum((PLM_diff.*PLM_visits).^2))/mysum_)^0.5;      %RMSdiff computed only over the visited points of the fine grid - weighted by number of visits to each point in that grid
    diff_vector = [diff_vector;diff];
    
    % This is just for the plotting routines
    for iN=1:nval_NN
        for iB=1:nval_BB
            if PLM_visits(iB,iN)==0
                PLM_visits(iB,iN)=NaN;   % substitute zeros for NaN so it looks better in the graphs
            end
        end
    end
    
    % Extra: Just to check that the NN works better thatn the regression, let's run a linear regression and calculate its RMSE and R^2
    
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
    surf(NN_grid,BB_grid,PLM_finegrid)
    title('Perceived Law of Motion','interpreter','Latex')
    xlabel('N','interpreter','Latex')
    ylabel('B','interpreter','Latex')
    axis([ Nmin Nmax Bmin Bmax -0.15 0.15])
    savefig(myfig,strcat('fig',num2str(it_PLM),'.fig'));
    close all
    
    for iN=1:nval_NN
        for iB=1:nval_BB
           if and(isnan(PLM_visits(iB,iN)),isnan(PLM_visits_(iB,iN)))
               PLM_diff(iB,iN)=NaN;                                   % put NaN on non-visited points, so they don't appear in the graph
           end
        end
    end
    
    myfig=figure(1);
    surf(NN_grid,BB_grid,PLM_diff);
    saveas(myfig,strcat('figdif',num2str(it_PLM),'.fig'))
    close all
    
    wePLM = wePLM*wePLM1+wePLM2;
    
    if and(diff<critPLM,it_PLM>5)  % stop if diff meets criterion, but not if we've ran surprisingly few iterations

        PLM = PLM_2;
        PLM_finegrid = PLM_finegrid_2;

        myfig=figure(1);
        surf(NN_grid,BB_grid,PLM_finegrid)
        title('Perceived Law of Motion','interpreter','Latex')
        xlabel('N','interpreter','Latex')
        ylabel('B','interpreter','Latex')
        axis([ Nmin Nmax Bmin Bmax -0.15 0.15])
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

