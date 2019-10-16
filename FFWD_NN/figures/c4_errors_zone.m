% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all



%%
% PLM forecast errors

eval(load_LR_errors);
Y_errors_LR = Y-Y_fit;

load '../z_FinalWorkspace.mat' Y Y_fit;
Y_errors_NN = Y-Y_fit;



myfig=figure(20);
set(myfig, 'Position', [0 0 800 400])

[y_ax,x_ax] = hist(Y_errors_LR,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'-.','Color',[1,0.1,0.1],'Linewidth',2)
hold on
xlim([-0.005 0.005])
title('Forecast errors', 'interpreter','latex','FontSize',14);
ax = gca;
ax.YAxis.Exponent = 0;

[y_ax,x_ax] = hist(Y_errors_NN,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'-','Color',[0.1,0.3,1],'Linewidth',2)
xlim([-0.005 0.005])
ax = gca;
ax.YAxis.Exponent = 0;
legend({'Linear','Neural network'}, 'interpreter','latex','FontSize',12)

grid

print -dpdf h20_PLMerrors
savefig(myfig,'h20_PLMerrors.fig');

title(' ', 'interpreter','latex','FontSize',14);

print -dpdf g20_PLMerrors




%%
% distribution of visited points

distrib_B = zeros(size(PLM_visits,1),1);
distrib_N = zeros(size(PLM_visits,2),1);
K_surf    = zeros(size(PLM_visits))*NaN;

PLM_visits_sum=0;
for iB=1:size(PLM_visits,1)
    for iN=1:size(PLM_visits,2)
        if not(isnan(PLM_visits(iB,iN)))
            PLM_visits_sum=PLM_visits_sum+PLM_visits(iB,iN);
            distrib_N(iN)=distrib_N(iN)+PLM_visits(iB,iN);
            distrib_B(iB)=distrib_B(iB)+PLM_visits(iB,iN);
            K_surf(iB,iN)=BB_grid(iB)+NN_grid(iN);
        end
    end
end

PLM_visits_ = PLM_visits ./ (PLM_visits_sum*dBB*dNN);
distrib_B_  = distrib_B  ./ (sum(distrib_B(:))*dBB);
distrib_N_  = distrib_N  ./ (sum(distrib_N(:))*dNN);

% 2D graph

myfig=figure(34);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,[1:4])
h=pcolor(squeeze(BB_grid),squeeze(NN_grid),PLM_visits');
set(h, 'EdgeColor', 'none');
hold on
plot(B_llsss,N_llsss,'ro')
plot(B_hlsss ,N_hlsss,'ro')
plot(B3 ,N3,'ro')
text(B_llsss,N_llsss,'"Low leverage SSS" $\rightarrow$', 'interpreter','latex','HorizontalAlignment','right','FontSize',12)
text(B_hlsss ,N_hlsss,'"Baseline SSS" $\rightarrow$', 'interpreter','latex','HorizontalAlignment','right','FontSize',12)
text(B3,N3,'Arbitrary high-leverage point $\rightarrow$', 'interpreter','latex','HorizontalAlignment','right','FontSize',12)
title('(a) Ergodic distribution $f(B,N)$', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
grid

subplot(3,2,5)
plot(BB_grid,distrib_B_,'-b','LineWidth',2)
title('(b) Marginal distribution of debt ($B$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
grid

subplot(3,2,6)
plot(NN_grid,distrib_N_,'-b','LineWidth',2)
title('(c) Marginal distribution of equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Nmin Nmax])
grid

print -dpdf h34_zone_2D
savefig(myfig,'h34_zone_2D.fig');

% 3D graph

myfig=figure(35);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,[1:4])
mesh(squeeze(NN_grid),squeeze(BB_grid),PLM_visits_)
title('(a) Ergodic density $f(B,N)$', 'interpreter','latex','FontSize',14);
xlabel('equity ($N$)', 'interpreter','latex','FontSize',14);
ylabel('debt ($B$)', 'interpreter','latex','FontSize',14);
xlim([Nmin Nmax])
ylim([Bmin Bmax])

subplot(3,2,5)
plot(BB_grid,distrib_B_,'-k','LineWidth',2)
hold on
plot(BB_grid,(BB_grid>B_llsss)*100-50,'--','Color',[1,0.1,0.1],'Linewidth',1)
plot(BB_grid,(BB_grid>B_hlsss)*100-50, '-','Color',[0,0.5,0.1],'Linewidth',1)
plot(BB_grid,(BB_grid>B_ss)*100-50   ,'-.','Color',[0.1,0.3,1],'Linewidth',1)
title('(b) Marginal distribution of debt ($B$)', 'interpreter','latex','FontSize',14);
legend({'Marginal distribution','LL-SSS','HL-SSS','DSS'},'Location','southoutside', 'interpreter','latex','FontSize',12)
xlim([Bmin Bmax])
ylim([0 3])
grid

subplot(3,2,6)
plot(NN_grid,distrib_N_,'-k','LineWidth',2)
hold on
plot(NN_grid,(NN_grid>N_llsss)*100-50,'--','Color',[1,0.1,0.1],'Linewidth',1)
plot(NN_grid,(NN_grid>N_hlsss)*100-50, '-','Color',[0,0.5,0.1],'Linewidth',1)
plot(NN_grid,(NN_grid>N_ss)*100-50   ,'-.','Color',[0.1,0.3,1],'Linewidth',1)
title('(c) Marginal distribution of equity ($N$)', 'interpreter','latex','FontSize',14);
legend({'Marginal distribution','LL-SSS','HL-SSS','DSS'},'Location','southoutside', 'interpreter','latex','FontSize',12)
xlim([Nmin Nmax])
ylim([0 1.5])
grid

print -dpdf h35_zone_3D
savefig(myfig,'h35_zone_3D.fig');

print -dpdf g35_zone_3D

