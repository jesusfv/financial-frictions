%MAIN PROGRAM
% Galo Nuno (2016)
clear all; clc;

%PARAMETERS
alpha         = 0.35;                       % Capital share   
delta         = 0.1;                        % Depreciation
gamma         = 2;                          % CRRA utility with parameter s
rho           = 0.05;                       % discount rate
la1           = 0.5;                        % transition prob
la2           = 1/3;
z1            = 0.93;                       % labor productivity
z2            = 1+ la2/la1 *(1- z1);
amin          = 0;                          % borrowing constraint
amax          = 30;                         % max value of assets
                       
I        = 1001;                            % numer of points between amin and amax, please make it end in 1
maxit    = 100;                             % max number of iterations
crit     = 10^(-6);                         % convergence ctrit 
Delta    = 1000;                            % Delta in HJB
maxitK   = 400;                             % max number of iterations
critK    = 10^(-3);                         % convergence ctrit
weK      = 0.001;                           % Weigth in the relaxation algorithm

parameters.alpha  = alpha;
parameters.delta  = delta;
parameters.gamma  = gamma;
parameters.rho    = rho;
parameters.z1     = z1;
parameters.z2     = z2;
parameters.la1    = la1;
parameters.la2    = la2;
parameters.I      = I;
parameters.amin   = amin;
parameters.amax   = amax;
parameters.maxit  = maxit;
parameters.crit   = crit;
parameters.Delta  = Delta;
parameters.maxitK  = maxitK;
parameters.critK   = critK;
parameters.weK     = weK;

%---------------------------------------------------
% STEADY-STATE
K_guess = 3.69;
K       = b2_equilibrium(parameters, K_guess);

results      = b3_HJB_stationary(parameters,K);
distribution = b4_KFE_stationary(parameters,results);
 
%OPTIONS
c1_Moments
c2_plotting

K_ss=S;
g_ss=distribution;
A_ss=results.A;
save 'ss_results.mat' A_ss K_ss g_ss -mat
