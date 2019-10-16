% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

% load parameters

run ('../b1_parameters')
load '../steadystate/ss_results'

% solve model using linear regression for PLM

tic
b2_Klm
disp('model solution using LR: completed')
toc

b9_plot

close all

% save 'z_FinalWorkspace_gsim.mat' g_big -mat -v7.3
clear g_big
save 'z_FinalWorkspace.mat' -mat -v7.3

