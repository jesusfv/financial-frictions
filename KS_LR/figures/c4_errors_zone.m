% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all



%%
% PLM forecast errors

Y_errors_LR = Y-Y_fit;



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


grid

print -dpdf h20_PLMerrors
savefig(myfig,'h20_PLMerrors.fig');

title(' ', 'interpreter','latex','FontSize',14);

print -dpdf g20_PLMerrors




%%
% distribution of visited points

distrib_K = zeros(size(PLM_visits,1),1);
distrib_Z = zeros(size(PLM_visits,2),1);

PLM_visits_sum=0;
for iK=1:size(PLM_visits,1)
    for iZ=1:size(PLM_visits,2)
        if not(isnan(PLM_visits(iK,iZ)))
            PLM_visits_sum=PLM_visits_sum+PLM_visits(iK,iZ);
            distrib_Z(iZ)=distrib_Z(iZ)+PLM_visits(iK,iZ);
            distrib_K(iK)=distrib_K(iK)+PLM_visits(iK,iZ);
        end
    end
end

PLM_visits_ = PLM_visits ./ (PLM_visits_sum*dKK*dZ);
distrib_K_  = distrib_K  ./ (sum(distrib_K(:))*dKK);
distrib_Z_  = distrib_Z  ./ (sum(distrib_Z(:))*dZ);

% 2D graph

myfig=figure(34);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,[1:4])
h=pcolor(squeeze(KK_grid),squeeze(Z_grid),PLM_visits');
set(h, 'EdgeColor', 'none');
hold on
plot(K_sss,Z_sss,'ro')
title('(a) Ergodic distribution $f(K,Z)$', 'interpreter','latex','FontSize',14);
xlabel('capital ($K$)', 'interpreter','latex','FontSize',14);
ylabel('productivity ($Z$)', 'interpreter','latex','FontSize',14);
xlim([Kmin Kmax])
ylim([Zmin Zmax])
grid

subplot(3,2,5)
plot(KK_grid,distrib_K_,'-b','LineWidth',2)
title('(b) Marginal distribution of capital ($K$)', 'interpreter','latex','FontSize',14);
xlim([Kmin Kmax])
grid

subplot(3,2,6)
plot(Z_grid,distrib_Z_,'-b','LineWidth',2)
title('(c) Marginal distribution of productivity ($Z$)', 'interpreter','latex','FontSize',14);
xlim([Zmin Zmax])
grid

print -dpdf h34_zone_2D
savefig(myfig,'h34_zone_2D.fig');

% 3D graph

myfig=figure(35);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,[1:4])
mesh(squeeze(Z_grid),squeeze(KK_grid),PLM_visits_)
title('(a) Ergodic density $f(K,Z)$', 'interpreter','latex','FontSize',14);
xlabel('capital ($K$)', 'interpreter','latex','FontSize',14);
ylabel('productivity ($Z$)', 'interpreter','latex','FontSize',14);
xlim([Zmin Zmax])
ylim([Kmin Kmax])

subplot(3,2,5)
plot(KK_grid,distrib_K_,'-k','LineWidth',2)
hold on
plot(KK_grid,(KK_grid>K_sss)*10000-5000,'--','Color',[1,0.1,0.1],'Linewidth',1)
plot(KK_grid,(KK_grid>K_ss)*10000-5000 ,'-.','Color',[0.1,0.3,1],'Linewidth',1)
title('(b) Marginal distribution of capial ($K$)', 'interpreter','latex','FontSize',14);
legend({'Marginal distribution','SSS','DSS'},'Location','southoutside', 'interpreter','latex','FontSize',12)
xlim([Kmin Kmax])
ylim([0 15])
grid

subplot(3,2,6)
plot(Z_grid,distrib_Z_,'-k','LineWidth',2)
hold on
plot(Z_grid,(Z_grid>Z_sss)*10000-5000,'--','Color',[1,0.1,0.1],'Linewidth',1)
plot(Z_grid,(Z_grid>Zmean)*10000-5000,'-.','Color',[0.1,0.3,1],'Linewidth',1)
title('(c) Marginal distribution of productivity ($Z$)', 'interpreter','latex','FontSize',14);
legend({'Marginal distribution','SSS','DSS'},'Location','southoutside', 'interpreter','latex','FontSize',12)
xlim([Zmin Zmax])
ylim([0 50])
grid

print -dpdf h35_zone_3D
savefig(myfig,'h35_zone_3D.fig');

print -dpdf g35_zone_3D

