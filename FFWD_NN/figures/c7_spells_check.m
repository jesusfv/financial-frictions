% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

% assign each simulated point to a basin

basin_MNR = mnrfit([basin_B basin_N],basin_ll+1);
Bsim = basin_Bsim_big(:);
Nsim = basin_Nsim_big(:);
basin_MNRfit = [ones(size(Bsim)) Bsim Nsim] * basin_MNR;

sim_zone = (basin_MNRfit >0); % one if the point belongs to the hlsss basin

clear basin_MNR basn_MNRfit

% plot results



Bsim_hlsss = Bsim .* (1./sim_zone);
Nsim_hlsss = Nsim .* (1./sim_zone);
Bsim_llsss = Bsim .* (1./(1-sim_zone));
Nsim_llsss = Nsim .* (1./(1-sim_zone));


myfig=figure(93);
set(myfig, 'Position', [0 0 800 800])
plot(Bsim_hlsss,Nsim_hlsss,'.','Color',[0.1,0.3,1],'linewidth',1);
hold on;
plot(Bsim_llsss,Nsim_llsss,'.','Color',[1,0.1,0.1],'linewidth',1);
title('Check that the multinomial logistic regression is doing a good job', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
