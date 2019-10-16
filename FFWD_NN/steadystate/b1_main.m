% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

clear all; clc;

run('../b1_parameters')


I        = nval_a;                          % number of points in amin-to-amax range (individual savings)
maxit    = 100;                             % max number of iterations
crit     = 10^(-6);                         % convergence crit 
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
parameters.maxitK = maxitK;
parameters.critK  = critK;
parameters.weK    = weK;

%---------------------------------------------------
% STEADY-STATE

K            = ((rhohat+delta)/alpha)^(1/(alpha-1)); 
results      = b3_HJB_stationary(parameters,K);
distribution = b4_KFE_stationary(parameters,results);
 
c1_Moments
c2_plotting

B_ss = S;
g_ss = distribution;
A_ss = results.A;
K_ss = K;
N_ss = K_ss - B_ss;
save 'ss_results.mat' A_ss B_ss g_ss N_ss K_ss -mat
