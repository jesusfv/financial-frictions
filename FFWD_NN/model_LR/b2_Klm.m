% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

% VARIABLES:
% the main ones will be 4-dimmensional matrices
% of size (nval_a, nval_z, nval_B, nval_N), with subscripts (ia, iz, iB, iN)
% there will be lots of repeated information in them (many of these matrices
% will be flat in all dimensions except one), but I can afford it
% (with 51x2x61x41, each of these matrices takes slightly over 2 MB of memory)

a_grid   = linspace(amin,amax,nval_a)';           % 1D - assets
z_grid   = [z1,z2];                               % 1D - productivity
B_grid   = linspace(Bmin,Bmax,nval_B)';           % 1D - capital
N_grid   = linspace(Nmin,Nmax,nval_N)';           % 1D - TFP

BB_grid  = linspace(Bmin,Bmax,nval_BB)';          % finer grid, used only for determining visited range and convergence criteria
NN_grid  = linspace(Nmin,Nmax,nval_NN)';          % finer grid, used only for determining visited range and convergence criteria

a      = zeros(nval_a, nval_z, nval_B, nval_N);
z      = zeros(nval_a, nval_z, nval_B, nval_N);
B      = zeros(nval_a, nval_z, nval_B, nval_N);
N      = zeros(nval_a, nval_z, nval_B, nval_N);

BB_grid_2D = zeros(nval_BB, nval_NN);
NN_grid_2D = zeros(nval_BB, nval_NN);

for iz=1:nval_z    % repmat would be faster, but this is clearer
    for iB=1:nval_B
        for iN=1:nval_N
            a(:,iz,iB,iN)=a_grid;
        end
    end
end

for ia=1:nval_a
    for iB=1:nval_B
        for iN=1:nval_N
            z(ia,:,iB,iN)=z_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iN=1:nval_N
            B(ia,iz,:,iN)=B_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iB=1:nval_B
            N(ia,iz,iB,:)=N_grid;
        end
    end
end

for iN=1:nval_NN
    BB_grid_2D(:,iN)=BB_grid;
end
for iB=1:nval_BB
    NN_grid_2D(iB,:)=NN_grid;
end

a2=squeeze(a(:,:,1,1));    % this one is 2D instead of 4D, we need it for a simpler KFE algorithm

% Interest rates and wages (4D matrices that don't depend on anything but parameters) - WE ARE ASSUMING L=1
r =  alpha * Zeta * ((B+N).^(alpha-1)) - delta - sigma2*((B+N)./N);
w = (1-alpha) * Zeta * (B+N).^alpha;

% INITIAL GUESS - by bringing it here, we make the HJB start from the previous V1 optimum every time it starts with a new PLM
V1      = ((w.*z + r.*a).^(1-gamma)-1)/(1-gamma)/rho;  % utility from always consuming today's income
% construct initial V guess in shape of V1_stacked
V1_stacked = zeros(nval_a*nval_z*nval_B*nval_N,1);

for iN=1:nval_N    % reshape could be faster, but this is clearer
    for iB=1:nval_B
        for iz=1:nval_z
            V1_stacked(1+(iz-1)*nval_a+(iB-1)*nval_a*nval_z+(iN-1)*nval_a*nval_z*nval_B:nval_a+(iz-1)*nval_a+(iB-1)*nval_a*nval_z+(iN-1)*nval_a*nval_z*nval_B,1)=V1(:,iz,iB,iN);
        end
    end
end


% random draws for KFE simulations (they are always the same)

rng(rngseed1);
e_shocks = randn(nval_sim,1);
sdte_shocks = sigma*(dt^0.5)*e_shocks;

% PLM iteration
PLM   = zeros(nval_a, nval_z, nval_B, nval_N);   % initial guess for tomorrow's B is today's B (zero change)
                                                 % (the PLM will always be flat on its first two dimmensions; but it's easier to have it as 4D)
PLM_2 = zeros(nval_a, nval_z, nval_B, nval_N);

PLM_finegrid   = zeros(nval_BB, nval_NN);        % finer grid, used only for determining visited range and PLM convergence criteria
PLM_finegrid_2 = zeros(nval_BB, nval_NN);
PLM_visits     = zeros(nval_BB, nval_NN).*NaN;   % Number of visits in each grid point in the current iteration

X0=    ones(multi_sim*used_sim,1);               % we will use this as constant in the just-for-reporting linear regression

format long;

for it_PLM=1:maxitPLM

    PLM_visits_=PLM_visits;   % Number of visits in each grid point in the previous iteration
    
    b3_HJB;
    b5_KFE;
    b7_PLM;

    t=0;
    PLM_visits=zeros(nval_BB,nval_NN);     % Let's count how many times we visited each point in the fine grid
    for it_sim=1:multi_sim
        t=t+delay_sim;
        for it_t=1:used_sim                % Looping over the actual simulation points, we will count "0.25 visits" for each of the four closest grid points to each simulated point
            t=t+1;
            PLM_visits(BBposD(t-1),NNposU(t-1))=PLM_visits(BBposD(t-1),NNposU(t-1))+0.25;
            PLM_visits(BBposU(t-1),NNposU(t-1))=PLM_visits(BBposU(t-1),NNposU(t-1))+0.25;
            PLM_visits(BBposD(t-1),NNposD(t-1))=PLM_visits(BBposD(t-1),NNposD(t-1))+0.25;
            PLM_visits(BBposU(t-1),NNposD(t-1))=PLM_visits(BBposU(t-1),NNposD(t-1))+0.25;
        end
    end

    for iN=1:nval_NN
        for iB=1:nval_BB
            if PLM_visits(iB,iN)==0
                PLM_visits(iB,iN)=NaN;
            end
        end
    end
    
    mysum=nval_NN*nval_BB;                                            % mysum will be the number of fine-grid points that have at least one visit; we start assuming all of them are visited
    PLM_diff=PLM_finegrid_2-PLM_finegrid;
    for iN=1:nval_NN                                
        for iB=1:nval_BB
           if and(isnan(PLM_visits(iB,iN)),isnan(PLM_visits_(iB,iN))) % if this grid point has not been visited in the current and previous set of simulations, this point is ignored
               PLM_diff(iB,iN)=0;
               mysum=mysum-1;                                         % and the number of visited points falls by one
           end
        end
    end
    diff = (sum(sum(PLM_diff.^2))/mysum)^0.5;      %RMSE computed only over the visited points of the fine grid
    
    % just to check that the NN works better, let's run a linear regression and calculate its RMSE and R^2
    
    X_LR = [X0 X1 X2];
    beta = (X_LR'*X_LR)^-1*X_LR'*Y;
    Y_LR_fit = X_LR*beta;
    Y_LR_RMSE = (mean((Y-Y_LR_fit).^2))^0.5;
    Y_LR_R_squared = 1- ( sum((Y-Y_LR_fit).^2) ) / ( sum((Y-mean(Y)).^2) );

    % report summary results from this iteration
    disp('   iter                Y_LR_RMSE           Y_LR_R_squared       diff_old:')
    disp([   it_PLM              Y_LR_RMSE           Y_LR_R_squared       diff]);
%   save 'PLM_iter.mat' Bsim Zsim PLM V1 c n -mat
%     eval(['save ''klm_iter_' num2str(n_klm) '.mat'' Ksim Zsim KK_lm V1 c n -mat'])
%     diff=0;
    PLM = (1-wePLM)*PLM + wePLM*PLM_2;
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
    
    if diff<critPLM

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
%%
