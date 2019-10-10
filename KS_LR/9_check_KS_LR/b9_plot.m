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
plot(Ksim)
savefig(myfig,'fig__1_Ksim.fig');

myfig=figure(2);
plot(Zsim)
savefig(myfig,'fig__2_Zsim.fig');

myfig=figure(4);
plot(Y)
hold on
plot(Y_fit);
title('Actual (blue) versus expected (red) growth rate of K','interpreter','Latex')
xlabel('Time (months)')
savefig(myfig,'fig__4_ACTUALvsPLM.fig');

myfig=figure(5);
plot(Y-Y_fit)
savefig(myfig,'fig__5_PLMerrors.fig');

myfig=figure(10);
histogram(Y-Y_fit)
xlim([-0.003 0.003])
savefig(myfig,'fig__5_PLMerrors_.fig');


%myfig=figure(6);
%plot(progress(:,minloc));
%savefig(myfig,'fig__6_NNtrain.fig');

%myfig=figure(7);
%semilogy(progress(:,minloc));
%savefig(myfig,'fig__7_NNtrainlog.fig');

myfig=figure(8);
surf(Z_grid,KK_grid,PLM_visits)
title('number of visits to each point in the grid')
xlabel('Z')
ylabel('K')
savefig(myfig,'fig__8_nvisits.fig');

myfig=figure(9);
surf(Z_grid,KK_grid,PLM_finegrid)
title('PLM (neural network)')
xlabel('Z')
ylabel('K')
savefig(myfig,'fig__9_PLM.fig');




