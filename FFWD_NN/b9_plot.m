% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

% plot some graphs that are useful for diagnosis

myfig=figure(11);
set(myfig, 'Position', [0 0 800 300])
subplot(1,2,1)
plot(diff_vector)
title('(a) Convergence', 'interpreter','latex','FontSize',12);
subplot(1,2,2)
plot(diff_vector)
ylim([0 0.002])
title('(b) Convergence (zoomed)', 'interpreter','latex','FontSize',12);
print -dpdf fig__0_convergence
print -depsc fig__0_convergence
savefig(myfig,'fig__0_convergence.fig');

myfig=figure(1);
plot(Bsim)
savefig(myfig,'fig__1_Bsim.fig');

myfig=figure(2);
plot(Nsim)
savefig(myfig,'fig__2_Nsim.fig');

myfig=figure(3);
plot(Bsim+Nsim)
savefig(myfig,'fig__3_Ksim.fig');

myfig=figure(4);
plot(Y)
hold on
plot(Y_fit);
title('Actual (blue) versus expected (red) growth rate of B','interpreter','Latex')
xlabel('Time (months)')
savefig(myfig,'fig__4_ACTUALvsPLM.fig');

myfig=figure(5);
plot(Y-Y_fit)
savefig(myfig,'fig__5_PLMerrors.fig');

myfig=figure(10);
histogram(Y-Y_fit)
xlim([-0.003 0.003])
savefig(myfig,'fig__5_PLMerrors_.fig');


myfig=figure(6);
plot(progress(:,minloc));
savefig(myfig,'fig__6_NNtrain.fig');

myfig=figure(7);
semilogy(progress(:,minloc));
savefig(myfig,'fig__7_NNtrainlog.fig');

myfig=figure(8);
surf(NN_grid,BB_grid,PLM_visits)
title('number of visits to each point in the grid')
xlabel('N')
ylabel('B')
savefig(myfig,'fig__8_nvisits.fig');

myfig=figure(9);
surf(NN_grid,BB_grid,PLM_finegrid)
title('PLM (neural network)')
xlabel('N')
ylabel('B')
savefig(myfig,'fig__9_PLM.fig');




% calculate standard deviation of simulated output growth rates

Ksim     = Bsim+Nsim   ;
Ysim     = Zeta * Ksim.^alpha   ;

YsimAg   = zeros(used_sim*multi_sim*dt,1);
t=0;
for it1=1:multi_sim
    t=t+delay_sim;
    for it2=1:used_sim*dt
        t=t+1/dt;
        YsimAg(t)=sum(Ysim(t-1/dt+1:t))/sum(Ysim(t-2/dt+1:t-1/dt))*100-100;
    end
end
YsimAg_desvest  = std(YsimAg)

YsimQg    = zeros(used_sim*multi_sim*dt*4,1);
t=0;
for it1=1:multi_sim
    t=t+delay_sim;
    for it2=1:used_sim*dt*4
        t=t+1/dt/4;
        YsimQg(t)=sum(Ysim(t-1/dt/4+1:t))/sum(Ysim(t-2/dt/4+1:t-1/dt/4))*100-100;
    end
end
YsimQg_desvest  = std(YsimQg)



