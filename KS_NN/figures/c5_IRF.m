% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

nval_IRF  = 50/dt;

shocksize = -2/sqrt(dt);


irf_Ksim  = zeros(nval_IRF,1);
irf_KposD = zeros(nval_IRF,1);
irf_KposU = zeros(nval_IRF,1);

irf_Zsim  = zeros(nval_IRF,1);
irf_ZposD = zeros(nval_IRF,1);
irf_ZposU = zeros(nval_IRF,1);
irf_wZ    = zeros(nval_IRF,1);

irf_gsim  = zeros(nval_a,nval_z,nval_IRF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% baseline scenario, at sss

e_irf = zeros(nval_IRF,1);
g1    = g_sss;

c5_IRF_sim

irf_sssA_K = irf_Ksim;
irf_sssA_Z = irf_Zsim;
irf_sssA_g = irf_gsim;

% alternative scenario with shock, at sss

e_irf(1) = shocksize;
g1       = g_sss;

c5_IRF_sim

irf_sssB_K = irf_Ksim;
irf_sssB_Z = irf_Zsim;
irf_sssB_g = irf_gsim;


% alternative scenario with big shock, at sss

e_irf(1) = shocksize*2;
g1       = g_sss;

c5_IRF_sim

irf_sssB2_K = irf_Ksim;
irf_sssB2_Z = irf_Zsim;
irf_sssB2_g = irf_gsim;


% alternative scenario with positive shock, at sss

e_irf(1) = -shocksize;
g1       = g_sss;

c5_IRF_sim

irf_sssBP_K = irf_Ksim;
irf_sssBP_Z = irf_Zsim;
irf_sssBP_g = irf_gsim;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE DIFFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IRF_sss_K    = irf_sssB_K    ./ irf_sssA_K    *100-100;
IRF_sss_Z    = irf_sssB_Z    -  irf_sssA_Z    ;
IRF_sss_g    = irf_sssB_g    -  irf_sssA_g    ;

IRF_sss2_K    = irf_sssB2_K  ./ irf_sssA_K    *100-100;
IRF_sss2_Z    = irf_sssB2_Z  -  irf_sssA_Z    ;
IRF_sss2_g    = irf_sssB2_g  -  irf_sssA_g    ;

IRF_sssP_K    = irf_sssBP_K  ./ irf_sssA_K    *100-100;
IRF_sssP_Z    = irf_sssBP_Z  -  irf_sssA_Z    ;
IRF_sssP_g    = irf_sssBP_g  -  irf_sssA_g    ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% NORMALIZE DIFFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IRF_sss2_K    = IRF_sss2_K/2;
IRF_sss2_Z    = IRF_sss2_Z/2;
IRF_sss2_g    = IRF_sss2_g/2;

IRF_sssP_K    = -IRF_sssP_K;
IRF_sssP_Z    = -IRF_sssP_Z;
IRF_sssP_g    = -IRF_sssP_g;



%%
%%%%%%%%%%%%%%%
% plot graphs %
%%%%%%%%%%%%%%%

%%
% basic IRF


myfig=figure(50);
set(myfig, 'Position', [0 0 800 500])

subplot(1,2,1);
plot([1:nval_IRF]*dt, IRF_sss_K,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_sss2_K,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_sssP_K,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(a) Capital, $K$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(1,2,2);
plot([1:nval_IRF]*dt, IRF_sss_Z,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_sss2_Z,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_sssP_Z,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(b) Productivity, $Z$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

legend({'negative shock','big negative shock (rescaled)','positive shock (rescaled)'},'Location','best', 'interpreter','latex','FontSize',10)
% legend boxoff

print -dpdf h50_IRF
savefig(myfig,'h50_IRF.fig');

print -dpdf g50_IRF


%%
% DIRF

ylimL = floor(min(min((squeeze(IRF_sss_g(1:100,1,1:1/dt:40/dt)+IRF_sss_g(1:100,2,1:1/dt:40/dt)))*100)))/100;
ylimH = ceil (max(max((squeeze(IRF_sss_g(1:100,1,1:1/dt:40/dt)+IRF_sss_g(1:100,2,1:1/dt:40/dt)))*100)))/100;


myfig=figure(57);
set(myfig, 'Position', [0 0 600 600])

mesh([1:1/dt:40/dt]*dt,a_grid(1:100),squeeze(IRF_sss_g(1:100,1,1:1/dt:40/dt)+IRF_sss_g(1:100,2,1:1/dt:40/dt)));
title('DIRF', 'interpreter','latex','FontSize',12);
axis([0 25 a_grid(1) a_grid(50) ylimL ylimH])
view([110, 10])
xlabel('years','interpreter','latex','FontSize',12)
ylabel('assets ($a$)','interpreter','latex','FontSize',12)

print -dpdf h57_DIRF
savefig(myfig,'h57_DIRF.fig');

print -dpdf g57_DIRF


%% distributions, value function V, and effect of the shock on value function V

% V at sss, before the shock

KposD =floor((irf_sssA_K(1)-Kmin)/dK)+1;
KposU = ceil((irf_sssA_K(1)-Kmin)/dK)+1;
wK=(K_grid(KposU)-irf_sssA_K(1))/dK;  % weight of posD

ZposD =floor((irf_sssA_Z(1)-Zmin)/dZ)+1;
ZposU = ceil((irf_sssA_Z(1)-Zmin)/dZ)+1;
wZ=(Z_grid(ZposU)-irf_sssA_Z(1))/dZ;  % weight of posD

irf_sss_V1 = squeeze( wK*wZ*V1(:,:,KposD,ZposD) + (1-wK)*wZ*V1(:,:,KposU,ZposD) + wK*(1-wZ)*V1(:,:,KposD,ZposU) + (1-wK)*(1-wZ)*V1(:,:,KposU,ZposU) );

% V at sss, after the shock

KposD =floor((irf_sssB_K(2)-Kmin)/dK)+1;
KposU = ceil((irf_sssB_K(2)-Kmin)/dK)+1;
wK=(K_grid(KposU)-irf_sssB_K(1))/dK;  % weight of posD

ZposD =floor((irf_sssB_Z(2)-Zmin)/dZ)+1;
ZposU = ceil((irf_sssB_Z(2)-Zmin)/dZ)+1;
wZ=(Z_grid(ZposU)-irf_sssB_Z(1))/dZ;  % weight of posD

irf_sss_V2 = squeeze( wK*wZ*V1(:,:,KposD,ZposD) + (1-wK)*wZ*V1(:,:,KposU,ZposD) + wK*(1-wZ)*V1(:,:,KposD,ZposU) + (1-wK)*(1-wZ)*V1(:,:,KposU,ZposU) );




% now plot it

myfig=figure(58);
set(myfig, 'Position', [0 0 830 1000])

subplot(3,2,1);
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_sss(2:ceil(size(g_ss,1)/2),1)),'-','Color',[0,0.5,0.1],'Linewidth',2)
hold on
plot(a_grid(1),squeeze(g_sss(1,1)),'o','Color',[0,0.5,0.1],'Linewidth',2)
title('(a) Distribution ($g$) of low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,2);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_sss(1:ceil(size(g_ss,1)/2),2)),'-','Color',[0,0.5,0.1],'Linewidth',2)
title('(b) Distribution ($g$) of high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid


subplot(3,2,3);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_sss_V1(1:ceil(size(g_ss,1)/2),1)),'-','Color',[0,0.5,0.1],'Linewidth',2)
title('(c) Value function ($V$) of low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,4);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_sss_V1(1:ceil(size(g_ss,1)/2),2)),'-','Color',[0,0.5,0.1],'Linewidth',2)
title('(d) Value function ($V$) of high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid


subplot(3,2,5);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_sss_V2(1:ceil(size(g_ss,1)/2),1)-irf_sss_V1(1:ceil(size(g_ss,1)/2),1)),'-','Color',[0,0.5,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),zeros(ceil(size(g_ss,1)/2),1), '-k','Linewidth',1)
title('(e) Effect of shock on $V$, low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,6);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_sss_V2(1:ceil(size(g_ss,1)/2),2)-irf_sss_V1(1:ceil(size(g_ss,1)/2),2)),'-','Color',[0,0.5,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),zeros(ceil(size(g_ss,1)/2),1), '-k','Linewidth',1)
title('(f) Effect of shock on $V$, high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

print -dpdf h58_g_V_DV
savefig(myfig,'h58_g_V_DV.fig');

print -dpdf g58_g_V_DV

