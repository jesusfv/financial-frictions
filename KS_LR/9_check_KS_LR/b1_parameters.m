% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha         = 0.35;                       % Capital share   
delta         = 0.1;                        % Depreciation
gamma         = 2;                          % CRRA utility with parameter s
rho           = 0.05;                       % discount rate
Zmean         = 0.00;                       % mean for TFP
theta         = 0.50;                       % AR for TFP
sigma         = 0.01;                       % sigma for TFP
la1           = 0.5;                        % transition prob
la2           = 1/3;
la            = [la1,la2];                  % vector of transition prob
z1            = 0.93;                       % labor productivity
z2            = 1+ la2/la1 *(1- z1);
amin          = 0;                          % borrowing constraint
amax          = 30;                         % max value of assets
Kmin          = 3.5;                        % relevant range for aggregate K
Kmax          = 3.9;
Zmin          = -0.04;                      % relevant range for aggregate Z
Zmax          =  0.04;
                       
nval_a        = 1001;                       % number of points between amin and amax, please make it end in 1
nval_z        = 2;                          % number of options for z (the idiosincratic shock)
nval_K        = 4;                          % number of points between Kmin and Kmax
nmul_K        = 5;                          % number of interpolated points to add between every contiguous pair of K values
nval_KK       = (nval_K-1)*nmul_K+1;        % number of points between Kmin and Kmax, in the interpolated version of A
nval_Z        = 41;                         % number of points between Zmin and Zmax, please make it end in 1
nval_Zsim     = 12000;                      % number of periods in simulation
delay_Zsim    = 1200;                       % number of initial periods in simulation that will not be used

dt       = 1/12;                            % size of t jump
da  = (amax-amin)/(nval_a-1);               % size of a jump
dK  = (Kmax-Kmin)/(nval_K-1);               % size of K jump
dKK = (Kmax-Kmin)/(nval_KK-1);              % size of K jump after interpolation in K
dZ  = (Zmax-Zmin)/(nval_Z-1);               % size of Z jump


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



%---------------------------------------------------
% LOAD STEADY-STATE
load b0_ss_results



% VARIABLES:
% the main ones will be 4-dimmensional matrices
% of size (nval_a, nval_z, nval_K, nval_Z), with subscripts (ia, iz, iK, kZ)
% there will be lots of repeated information in them (many of these matrices
% will be flat in all dimensions except one), but I can afford it
% (with 51x2x61x41, each of these matrices takes slightly over 2 MB of memory)

a_grid   = linspace(amin,amax,nval_a)';           % 1D - assets
z_grid   = [z1,z2];                               % 1D - productivity
K_grid   = linspace(Kmin,Kmax,nval_K)';           % 1D - capital
KK_grid  = linspace(Kmin,Kmax,nval_KK)';          % 1D - capital, with interpolated values
Z_grid   = linspace(Zmin,Zmax,nval_Z)';           % 1D - TFP

a      = zeros(nval_a, nval_z, nval_K, nval_Z);
z      = zeros(nval_a, nval_z, nval_K, nval_Z);
K      = zeros(nval_a, nval_z, nval_K, nval_Z);
KK     = zeros(nval_a, nval_z, nval_KK,nval_Z);
Z      = zeros(nval_a, nval_z, nval_K, nval_Z);
ZZ     = zeros(nval_a, nval_z, nval_KK,nval_Z);

for iz=1:nval_z    % repmat would be faster, but this is clearer
    for iK=1:nval_K
        for iZ=1:nval_Z
            a(:,iz,iK,iZ)=a_grid;
        end
    end
end

for ia=1:nval_a
    for iK=1:nval_K
        for iZ=1:nval_Z
            z(ia,:,iK,iZ)=z_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iZ=1:nval_Z
            K(ia,iz,:,iZ)=K_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iZ=1:nval_Z
            KK(ia,iz,:,iZ)=KK_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iK=1:nval_K
            Z(ia,iz,iK,:)=Z_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iK=1:nval_KK
            ZZ(ia,iz,iK,:)=Z_grid;
        end
    end
end

KK_grid_2D = zeros(nval_KK, nval_Z);
Z_grid_2D = zeros(nval_KK, nval_Z);
for iZ=1:nval_Z
    KK_grid_2D(:,iZ)=KK_grid;
end
for iK=1:nval_KK
    Z_grid_2D(iK,:)=Z_grid;
end

a2=squeeze(a(:,:,1,1));    % this one is 2D instead of 4D, we need it for a simpler KFE algorithm

% Interest rates and wages (4D matrices that don't depend on anything but parameters) - WE ARE ASSUMING L=1
r =  alpha * K.^(alpha-1).*exp(Z) - delta;
w = (1-alpha) * K.^alpha .*exp(Z);
