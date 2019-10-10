% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha         = 0.35;                       % Capital share   
Zeta          = 1;                          % Aggregate productivity
delta         = 0.1;                        % Depreciation
gamma         = 2;                          % CRRA utility with parameter s
rho           = 0.05;                       % discount rate
rhohat        = 0.04971;                    % discount rate for banks
la1           = 0.986;                      % transition probability
la2           = 0.052;
la            = [la1,la2];                  % vector of transition prob

amin          = 0;                          % borrowing constraint
amax          = 20;                         % max value of individual savings
z1            = 0.72;                       % labor productivity
z2            = 1 + la2/la1 * (1-z1);

Bmin          = 0.7;                        % relevant range for aggregate savings
Bmax          = 2.7;
Nmin          = 1.2;                        % relevant range for aggregate equity
Nmax          = 3.2;
                       
nval_a        = 501;                        % number of points in amin-to-amax range (individual savings)
nval_z        = 2;                          % number of options for z (the idiosincratic shock)
nval_B        = 4;                          % number of points in Bmin-to-Bmax range (aggregate savings), on the coarse grid used for the HJB
nval_N        = 51;                         % number of points in Nmin-to-Nmax range (aggregate equity) , on the coarse grid used for the HJB

nval_BB       = 101;                        % finer grid, used for training the NN, for determining visited range and for the convergence criteria
nval_NN       = 101;                        % finer grid, used for training the NN, for determining visited range and for the convergence criteria

dt            = 1/12;                       % size of t jump
da  = (amax-amin)/(nval_a-1);               % size of a jump
dB  = (Bmax-Bmin)/(nval_B-1);               % size of B jump on the coarse grid
dN  = (Nmax-Nmin)/(nval_N-1);               % size of N jump on the coarse grid

dBB = (Bmax-Bmin)/(nval_BB-1);              % size of B jump on the fine grid
dNN = (Nmax-Nmin)/(nval_NN-1);              % size of B jump on the fine grid

sigma         = 0.0140;                     % sigma for aggregate capital quality shock
sigma2        = sigma^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multi_sim    =   4    ;                     % number of simulation starts (i.e. one plus number of times we stop the simulation and start again from the ss)
delay_sim    =  500/dt;                     % number of initial periods in each simulation that will not be used (pre-heat)
used_sim     = 5000/dt;                     % number of periods in each simulation that will actually be used
each_sim     = delay_sim+used_sim;          % number of periods on each simulation
nval_sim     = multi_sim*each_sim;          % total number of periods in all simulations
rngseed1     = 123;                         % RNG seed for calculating the shocks for the simulation
                                            % there used to be a second seed for the mini-batch selection in the NN training, but we got rid of that by moving to batch gradient descent

maxitHJB = 100;                             % max number of iterations for the HJB
critHJB  = 10^(-6);                         % convergence crit for the HJB 
weHJB    = 0.5;                             % relaxation algorithm for HJB
Delta    = 1000;                            % Delta in HJB
maxitPLM = 200;                             % max number of iterations of the full algorithm
critPLM  = 0.00050;                         % convergence crit for determining that the PLM has converged
wePLM    = 0.3;                             % Initial weigth in the relaxation algorithm for PLM convergence
wePLM1   = 0.9;                             % reduction of the relaxation algorithm: wePLM = wePLM*wePLM1+wePLM2
wePLM2   = 0.005;                           % reduction of the relaxation algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural network parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
network_width  =    16;                     % Number of neurons in the hidden layer
mynoise        =     1;                     % Size of random initial NN parameters
lambda         =   0.1;                     % NN regularization parameter
NN_iters       = 10000;                     % Number of iterations to train the network
NN_starts      =    10;                     % Number of random restarts of the network training (only affects the first step, then it just starts from the previous NN)
learning_speed =  0.01;                     % Only affects the first step, then this becomes adaptative
reglimY        =     4;                     % NN input normalization
reglimX        =     4;                     % NN input normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% VARIABLES:
% the main ones will be 4-dimmensional matrices
% of size (nval_a, nval_z, nval_B, nval_N), with subscripts (ia, iz, iB, iN)
% there will be lots of repeated information in them (many of these matrices
% will be flat in all dimensions except one), but I can afford it
% (in terms of memory these matrices are not that big anyway)

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
