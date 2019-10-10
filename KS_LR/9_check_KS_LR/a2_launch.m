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

cd steadystate
a2_launch
cd ..

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% special stuff just for this folder
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

b1_parameters;


%% %%%%%%% %%
% solve model
%% %%%%%%% %%

tic
b2_Klm
disp('model solution: completed')
toc

% plot and save results

warning('off','all')

b9_plot

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


