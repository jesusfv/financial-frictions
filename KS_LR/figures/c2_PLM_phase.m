% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use LR to calculate the PLM on the fine grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1m=zeros(nval_KK*nval_Z,1);
X2m=zeros(nval_KK*nval_Z,1);

for iZ=1:nval_Z
    for iK=1:nval_KK
        X1m(iZ*nval_KK-nval_KK+iK,1)=KK_grid(iK);
        X2m(iZ*nval_KK-nval_KK+iK,1)=Z_grid(iZ);
    end
end


X0m = ones(size(X1m));
X1m=log(X1m);
X3m=X1m.*X2m;

X_LR = [X0m X1m X2m X3m];
Y_mat2 = X_LR*beta;
PLM_surf= reshape(Y_mat2,[nval_KK nval_Z]);

% visited part of that PLM
PLM_surf_v = PLM_surf;
for iZ=1:nval_Z
    for iK=1:nval_KK
       if isnan(PLM_visits(iK,iZ))
           PLM_surf_v(iK,iZ)=NaN;      % put NaN on non-visited points
       end
    end
end




%%
%%%%%%%%%%%%%
% find PLM0 %
%%%%%%%%%%%%%

PLM0=NaN(nval_KK,2);

for it1=1:nval_KK
    found=0;
    it3 = NaN;
    for it2=1:nval_Z
        if PLM_surf(it1,it2)<0
            it3=it2;
            found=1;
        end
    end
    if it3>=nval_Z
        found = 0;
    end
    if found == 1
        myweight=PLM_surf(it1,it3+1)/(PLM_surf(it1,it3+1)-PLM_surf(it1,it3));
        PLM0(it1,1)=KK_grid(it1);
        PLM0(it1,2)=Z_grid(it3)*myweight+Z_grid(it3+1)*(1-myweight);
    end
end


%%%%%%%%%%%%%%
% find ZPLM0 %
%%%%%%%%%%%%%%

ZPLM0=NaN(nval_KK,2);

for it1=1:nval_KK
    ZPLM0(it1,1)=KK_grid(it1);
    ZPLM0(it1,2)=Zmean;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all Stochastic Steady States %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find SSS points (it's just for graphs, so: using linear approximation around zero of PLM and NPLM surfaces)
% Later I'll use a very long no-shocks simulation to refine these points, before using them to run IRFs etc
PLM0dif=PLM0;
PLM0dif(:,2)=PLM0(:,2)-ZPLM0(:,2);
SSS_points=[];
for it1=2:size(PLM0,1)
    if (PLM0dif(it1,2)*PLM0dif(it1-1,2))<0
        myweight=PLM0dif(it1,2)/(PLM0dif(it1,2)-PLM0dif(it1-1,2)); % weight of the it1-1 point
        SSS_K = myweight*PLM0dif(it1-1,1)+(1-myweight)*PLM0dif(it1,1);
        SSS_Z1 = myweight* PLM0(it1-1,2)+(1-myweight)* PLM0(it1,2);
        SSS_Z2 = myweight*ZPLM0(it1-1,2)+(1-myweight)*ZPLM0(it1,2);
        SSS_points=[SSS_points ; [SSS_K (SSS_Z1+SSS_Z2)/2] ];
    end
end

K_sss = SSS_points(1,1);
Z_sss = SSS_points(1,2);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use LR to calculate the PLM on specific cuts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cut at K_sss, full range

X1m = K_sss*ones(size(Z_grid));
X2m = Z_grid;

X0m = ones(size(X1m));
X1m=log(X1m);
X3m=X1m.*X2m;

X_LR = [X0m X1m X2m X3m];
PLM_Ksss = X_LR*beta;



% cut at Z_sss, full range

X1m = KK_grid;
X2m = Zmean*ones(size(KK_grid));

X0m = ones(size(X1m));
X1m=log(X1m);
X3m=X1m.*X2m;

X_LR = [X0m X1m X2m X3m];
PLM_Zsss = X_LR*beta;



%%
%%%%%%%%%%%%%%%%%%%%%%%
% Now plot everything %
%%%%%%%%%%%%%%%%%%%%%%%

%%
% PLM

myfig=figure(10);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,1);
plot(Z_grid,PLM_Ksss,'-','Color',[1,0.1,0.1]);
title('(a) PLM when $K=K_sss$', 'interpreter','latex','FontSize',14);
xlabel('shock ($Z$)', 'interpreter','latex','FontSize',14);
ylabel('PLM', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,2);
plot(KK_grid,PLM_Zsss,'-','Color',[1,0.1,0.1]);
title('(a) PLM when $Z=Zmean$', 'interpreter','latex','FontSize',14);
xlabel('capital ($K$)', 'interpreter','latex','FontSize',14);
ylabel('PLM', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,[3:6]);
mesh(KK_grid,Z_grid,PLM_surf');
hold on
surf(KK_grid,Z_grid,PLM_surf_v');
plot(PLM0(:,1),PLM0(:,2),'-','Color',[1,0.1,0.1]);
title('(c) The perceived law of motion, PLM', 'interpreter','latex','FontSize',14);
xlabel('capital ($K$)', 'interpreter','latex','FontSize',14);
ylabel('shock ($Z$)', 'interpreter','latex','FontSize',14);
xlim([Kmin Kmax])
ylim([Zmin Zmax])
view([60, 15])

print -dpdf h10_PLM
savefig(myfig,'h10_PLM.fig');

print -dpdf g10_PLM



%%
% phase diagram

myfig=figure(16);
set(myfig, 'Position', [0 0 800 800])

plot(PLM0(:,1),PLM0(:,2),'-','Color',[0.1,0.3,1],'linewidth',2);
hold on;
plot(ZPLM0(:,1),ZPLM0(:,2),'--r','linewidth',2);
title('Phase diagram', 'interpreter','latex','FontSize',14);
xlabel('capital ($K$)', 'interpreter','latex','FontSize',14);
ylabel('shock ($Z$)', 'interpreter','latex','FontSize',14);
xlim([Kmin Kmax])
ylim([Zmin Zmax])


plot(SSS_points(:,1),SSS_points(:,2),' Ok','linewidth',2)
plot(K_ss,Zmean,' s','Color',[0,0.5,0.1],'linewidth',5,'markers',10)

legend({'PLM=0','Z=Zmean$','Stochastic steady state','Deterministic steady state'},'Location','southwest', 'interpreter','latex','FontSize',12)
grid

print -dpdf h30_phase
savefig(myfig,'h30_phase.fig');

title(' ');
print -dpdf g30_phase



