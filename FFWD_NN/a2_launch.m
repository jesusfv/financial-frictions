% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This script launches the routines in charge of solving the model, plotting the results and saving them

clear all
close all
clc

diary off;
delete('_log.txt');
diary('_log.txt');

% calculate steady state
disp(' ')
disp('calculate steady state')
disp(' ')

cd steadystate
a2_launch
close all
cd ..

% solve model using linear regression
disp(' ')
disp('solve model using linear regression')
disp(' ')

cd model_LR
a2_launch
close all
cd ..


% make sure important stuff is in memory

b1_parameters;
load './steadystate/ss_results';

load_PLM       = 'load  ''./model_LR/z_FinalWorkspace.mat'' PLM PLM_2 PLM_finegrid PLM_finegrid_2 PLM_visits;';
load_LR_errors = 'load ''../model_LR/z_FinalWorkspace.mat'' Y Y_fit;';


%% %%%%%%% %%
% solve model
%% %%%%%%% %%
disp(' ')
disp('solve model using neural network')
disp(' ')

tic
b2_Klm
disp('model solution using NN: completed')
toc

% plot and save results

warning('off','all')

b9_plot

save 'z_FinalWorkspace.mat' Y Y_fit Bsim Nsim -mat -v7.3

cd figures
a2_launch
close all
cd ..

warning('on','all')

close all

save 'z_FinalWorkspace_gsim.mat' g_big -mat -v7.3
clear g_big
save 'z_FinalWorkspace.mat' -mat -v7.3

diary off;


