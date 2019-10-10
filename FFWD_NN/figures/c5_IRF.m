% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

nval_IRF=100/dt;

shocksize = -2/sqrt(dt);


Birf     = zeros(nval_IRF,1);
BirfposU = zeros(nval_IRF,1);
BirfposD = zeros(nval_IRF,1);

Nirf     = zeros(nval_IRF,1);
NirfposU = zeros(nval_IRF,1);
NirfposD = zeros(nval_IRF,1);

cirf     = zeros(nval_a,nval_z,nval_IRF);
Cirf     = zeros(nval_IRF,1);
rirf     = zeros(nval_IRF,1);

girf     = zeros(nval_a,nval_z,nval_IRF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% baseline scenario, at hlsss

e_irf = zeros(nval_IRF,1);
g1=g_hlsss;
Nirf(1)=N_hlsss;

c5_IRF_sim

irf_hlsssA_B = Birf;
irf_hlsssA_N = Nirf;
irf_hlsssA_r = rirf;
irf_hlsssA_g = girf;
irf_hlsssA_c = cirf;
irf_hlsssA_C = Cirf;

% alternative scenario with shock, at hlsss

e_irf(1) = shocksize;
g1=g_hlsss;
Nirf(1)=N_hlsss;

c5_IRF_sim

irf_hlsssB_B = Birf;
irf_hlsssB_N = Nirf;
irf_hlsssB_r = rirf;
irf_hlsssB_g = girf;
irf_hlsssB_c = cirf;
irf_hlsssB_C = Cirf;


% baseline scenario, at llsss

e_irf = zeros(nval_IRF,1);
g1=g_llsss;
Nirf(1)=N_llsss;

c5_IRF_sim

irf_llsssA_B = Birf;
irf_llsssA_N = Nirf;
irf_llsssA_r = rirf;
irf_llsssA_g = girf;
irf_llsssA_c = cirf;
irf_llsssA_C = Cirf;

% alternative scenario with shock, at llsss

e_irf(1) = shocksize;
g1=g_llsss;
Nirf(1)=N_llsss;

c5_IRF_sim

irf_llsssB_B = Birf;
irf_llsssB_N = Nirf;
irf_llsssB_r = rirf;
irf_llsssB_g = girf;
irf_llsssB_c = cirf;
irf_llsssB_C = Cirf;

% baseline scenario at arbitrary point with high leverage

mydistance = sqrt((Bsim-B3).^2 + (Nsim-N3).^2);
[mydistance_,myindex]=min(mydistance);
clear mydistance mydistance_

e_irf = zeros(nval_IRF,1);
g1=squeeze(g_big(:,:,myindex));
Nirf(1)=Nsim(myindex);

c5_IRF_sim

irf_hlA_B = Birf;
irf_hlA_N = Nirf;
irf_hlA_r = rirf;
irf_hlA_g = girf;
irf_hlA_c = cirf;
irf_hlA_C = Cirf;

% alternative scenario with shock, at arbitrary point with high leverage

e_irf(1) = shocksize;
g1=squeeze(g_big(:,:,myindex));
Nirf(1)=Nsim(myindex);

c5_IRF_sim

irf_hlB_B = Birf;
irf_hlB_N = Nirf;
irf_hlB_r = rirf;
irf_hlB_g = girf;
irf_hlB_c = cirf;
irf_hlB_C = Cirf;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% CALCULATE SECONDARY VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

irf_hlsssA_K     = irf_hlsssA_B+irf_hlsssA_N   ;
irf_hlsssA_Y     = Zeta * irf_hlsssA_K.^alpha   ;
irf_hlsssA_w     = (1-alpha)*irf_hlsssA_Y   ;
irf_hlsssA_rk    = alpha*irf_hlsssA_Y./irf_hlsssA_K - delta   ;
irf_hlsssA_Chat  = rhohat*irf_hlsssA_N   ;
irf_hlsssA_KF    = [irf_hlsssA_K(2:end) ; irf_hlsssA_K(end)];

irf_hlsssB_K     = irf_hlsssB_B+irf_hlsssB_N   ;
irf_hlsssB_Y     = Zeta * irf_hlsssB_K.^alpha   ;
irf_hlsssB_w     = (1-alpha)*irf_hlsssB_Y   ;
irf_hlsssB_rk    = alpha*irf_hlsssB_Y./irf_hlsssB_K - delta   ;
irf_hlsssB_Chat  = rhohat*irf_hlsssB_N   ;
irf_hlsssB_KF    = [irf_hlsssB_K(2:end) ; irf_hlsssB_K(end)];

irf_llsssA_K     = irf_llsssA_B+irf_llsssA_N   ;
irf_llsssA_Y     = Zeta * irf_llsssA_K.^alpha   ;
irf_llsssA_w     = (1-alpha)*irf_llsssA_Y   ;
irf_llsssA_rk    = alpha*irf_llsssA_Y./irf_llsssA_K - delta   ;
irf_llsssA_Chat  = rhohat*irf_llsssA_N   ;
irf_llsssA_KF    = [irf_llsssA_K(2:end) ; irf_llsssA_K(end)];

irf_1lsssB_K     = irf_llsssB_B+irf_llsssB_N   ;
irf_llsssB_Y     = Zeta * irf_1lsssB_K.^alpha   ;
irf_llsssB_w     = (1-alpha)*irf_llsssB_Y   ;
irf_llsssB_rk    = alpha*irf_llsssB_Y./irf_1lsssB_K - delta   ;
irf_llsssB_Chat  = rhohat*irf_llsssB_N   ;
irf_llsssB_KF    = [irf_1lsssB_K(2:end) ; irf_1lsssB_K(end)];

irf_hlA_K     = irf_hlA_B+irf_hlA_N   ;
irf_hlA_Y     = Zeta * irf_hlA_K.^alpha   ;
irf_hlA_w     = (1-alpha)*irf_hlA_Y   ;
irf_hlA_rk    = alpha*irf_hlA_Y./irf_hlA_K - delta   ;
irf_hlA_Chat  = rhohat*irf_hlA_N   ;
irf_hlA_KF    = [irf_hlA_K(2:end) ; irf_hlA_K(end)];

irf_hlB_K     = irf_hlB_B+irf_hlB_N   ;
irf_hlB_Y     = Zeta * irf_hlB_K.^alpha   ;
irf_hlB_w     = (1-alpha)*irf_hlB_Y   ;
irf_hlB_rk    = alpha*irf_hlB_Y./irf_hlB_K - delta   ;
irf_hlB_Chat  = rhohat*irf_hlB_N   ;
irf_hlB_KF    = [irf_hlB_K(2:end) ; irf_hlB_K(end)];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE DIFFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IRF_hlsss_Y    = irf_hlsssB_Y    ./ irf_hlsssA_Y    *100-100;
IRF_hlsss_c    = irf_hlsssB_c    ./ irf_hlsssA_c    *100-100;
IRF_hlsss_C    = irf_hlsssB_C    ./ irf_hlsssA_C    *100-100;
IRF_hlsss_Chat = irf_hlsssB_Chat ./ irf_hlsssA_Chat *100-100;
IRF_hlsss_K    = irf_hlsssB_K    ./ irf_hlsssA_K    *100-100;
IRF_hlsss_B    = irf_hlsssB_B    ./ irf_hlsssA_B    *100-100;
IRF_hlsss_N    = irf_hlsssB_N    ./ irf_hlsssA_N    *100-100;
IRF_hlsss_w    = irf_hlsssB_w    ./ irf_hlsssA_w    *100-100;
IRF_hlsss_r    = irf_hlsssB_r    -  irf_hlsssA_r    ;
IRF_hlsss_rk   = irf_hlsssB_rk   -  irf_hlsssA_rk   ;
IRF_hlsss_g    = irf_hlsssB_g    -  irf_hlsssA_g    ;


IRF_llsss_Y    = irf_llsssB_Y    ./ irf_llsssA_Y    *100-100;
IRF_llsss_c    = irf_llsssB_c    ./ irf_llsssA_c    *100-100;
IRF_llsss_C    = irf_llsssB_C    ./ irf_llsssA_C    *100-100;
IRF_llsss_Chat = irf_llsssB_Chat ./ irf_llsssA_Chat *100-100;
IRF_llsss_K    = irf_1lsssB_K    ./ irf_llsssA_K    *100-100;
IRF_llsss_B    = irf_llsssB_B    ./ irf_llsssA_B    *100-100;
IRF_llsss_N    = irf_llsssB_N    ./ irf_llsssA_N    *100-100;
IRF_llsss_w    = irf_llsssB_w    ./ irf_llsssA_w    *100-100;
IRF_llsss_r    = irf_llsssB_r    -  irf_llsssA_r    ;
IRF_llsss_rk   = irf_llsssB_rk   -  irf_llsssA_rk   ;
IRF_llsss_g    = irf_llsssB_g    -  irf_llsssA_g    ;


IRF_hl_Y    = irf_hlB_Y    ./ irf_hlA_Y    *100-100;
IRF_hl_c    = irf_hlB_c    ./ irf_hlA_c    *100-100;
IRF_hl_C    = irf_hlB_C    ./ irf_hlA_C    *100-100;
IRF_hl_Chat = irf_hlB_Chat ./ irf_hlA_Chat *100-100;
IRF_hl_K    = irf_hlB_K    ./ irf_hlA_K    *100-100;
IRF_hl_B    = irf_hlB_B    ./ irf_hlA_B    *100-100;
IRF_hl_N    = irf_hlB_N    ./ irf_hlA_N    *100-100;
IRF_hl_w    = irf_hlB_w    ./ irf_hlA_w    *100-100;
IRF_hl_r    = irf_hlB_r    -  irf_hlA_r    ;
IRF_hl_rk   = irf_hlB_rk   -  irf_hlA_rk   ;
IRF_hl_g    = irf_hlB_g    -  irf_hlA_g    ;




%%
%%%%%%%%%%%%%%%
% plot graphs %
%%%%%%%%%%%%%%%

%%
% basic IRF


myfig=figure(50);
set(myfig, 'Position', [0 0 800 800])

subplot(3,3,1);
plot([1:nval_IRF]*dt, IRF_hlsss_Y,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_Y,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_Y,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(a) Output, $Y$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,2);
plot([1:nval_IRF]*dt, IRF_hlsss_C,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_C,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_C,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(b) HH consumption, $c$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,3);
plot([1:nval_IRF]*dt, IRF_hlsss_Chat,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_Chat,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_Chat,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(c) Expert consumption, $\hat{c}$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,4);
plot([1:nval_IRF]*dt, IRF_hlsss_K,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_K,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_K,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(d) Capital, $K$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,5);
plot([1:nval_IRF]*dt, IRF_hlsss_B,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_B,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_B,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(e) Debt, $B$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,6);
plot([1:nval_IRF]*dt, IRF_hlsss_N,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_N,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_N,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(f) Equity, $N$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,7);
plot([1:nval_IRF]*dt, IRF_hlsss_w,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_w,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_w,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(g) Wages, $w$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

legend({'HL-SSS','LL-SSS','HL point'},'Location','best', 'interpreter','latex','FontSize',10)
% legend boxoff

subplot(3,3,8);
plot([1:nval_IRF]*dt, IRF_hlsss_r,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_r,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_r,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(h) Interest rate, $r$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,9);
plot([1:nval_IRF]*dt, IRF_hlsss_rk,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_llsss_rk,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:nval_IRF]*dt, IRF_hl_rk,'--','Color',[0.1,0.3,1],'linewidth',2);
title('(i) Return to capital, $r_k$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

print -dpdf h50_IRF
savefig(myfig,'h50_IRF.fig');

print -dpdf g50_IRF

disp(' ')
disp('Present discounted value of shock to wages:')
disp(['LL-SSS: ' num2str(sum(exp(-rho*[0:nval_IRF-1]).*IRF_llsss_w(1:nval_IRF)'*dt))])
disp(['HL-SSS: ' num2str(sum(exp(-rho*[0:nval_IRF-1]).*IRF_hlsss_w(1:nval_IRF)'*dt))])
disp(['HL point: ' num2str(sum(exp(-rho*[0:nval_IRF-1]).*IRF_hl_w(1:nval_IRF)'*dt))])

IRF_llsss_w_025 = 0;
IRF_hlsss_w_025 = 0;
IRF_hl_w_025 = 0;
for it1 = 2:nval_IRF
    if IRF_llsss_w_025 == 0 & IRF_llsss_w(it1)>0.25*IRF_llsss_w(2)
        IRF_llsss_w_025 = it1*dt
    end
    if IRF_hlsss_w_025 == 0 & IRF_hlsss_w(it1)>0.25*IRF_hlsss_w(2)
        IRF_hlsss_w_025 = it1*dt
    end
    if IRF_hl_w_025 == 0 & IRF_hl_w(it1)>0.25*IRF_hl_w(2)
        IRF_hl_w_025 = it1*dt
    end
end
disp(' ')
disp('Years until less than 25% of the initial effect on wages is still alive:')
disp(['LL-SSS: ' num2str(IRF_llsss_w_025)])
disp(['HL-SSS: ' num2str(IRF_hlsss_w_025)])
disp(['HL point: ' num2str(IRF_hl_w_025)])


%%
% cIRF

ylimL1 = floor(min(min(squeeze(IRF_hlsss_c(1:100,2,1:20/dt))))*1)/1;
ylimH1 = ceil (max(min(squeeze(IRF_hlsss_c(1:100,2,1:20/dt))))*1)/1;

ylimL2 = floor(min(min(squeeze(IRF_llsss_c(1:100,2,1:20/dt))))*1)/1;
ylimH2 = ceil (max(min(squeeze(IRF_llsss_c(1:100,2,1:20/dt))))*1)/1;

ylimL = min([ylimL1 ylimL2]);
ylimH = max([ylimH1 ylimH2]);


myfig=figure(55);
set(myfig, 'Position', [0 0 800 500])

subplot(2,2,1);
plot(a_grid(1:100),squeeze(IRF_hlsss_c(1:100,2,2)),'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot(a_grid(1:100),squeeze(IRF_llsss_c(1:100,2,2)),'--','Color',[1,0.1,0.1],'linewidth',2);
title('(a) On impact', 'interpreter','latex','FontSize',12);
ylim([ylimL ylimH]);
ax = gca;
ax.YAxis.Exponent = 0;
grid on


subplot(2,2,2);
plot(a_grid(1:100),squeeze(IRF_hlsss_c(1:100,2,5/dt)),'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot(a_grid(1:100),squeeze(IRF_llsss_c(1:100,2,5/dt)),'--','Color',[1,0.1,0.1],'linewidth',2);
title('(b) After 5 years', 'interpreter','latex','FontSize',12);
ylim([ylimL ylimH]);
ax = gca;
ax.YAxis.Exponent = 0;
grid on


subplot(2,2,3);
plot(a_grid(1:100),squeeze(IRF_hlsss_c(1:100,2,10/dt)),'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot(a_grid(1:100),squeeze(IRF_llsss_c(1:100,2,10/dt)),'--','Color',[1,0.1,0.1],'linewidth',2);
title('(c) After 10 years', 'interpreter','latex','FontSize',12);
xlabel('assets ($a$)','interpreter','latex','FontSize',12)
ylim([ylimL ylimH]);
ax = gca;
ax.YAxis.Exponent = 0;
grid on


subplot(2,2,4);
plot(a_grid(1:100),squeeze(IRF_hlsss_c(1:100,2,20/dt)),'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot(a_grid(1:100),squeeze(IRF_llsss_c(1:100,2,20/dt)),'--','Color',[1,0.1,0.1],'linewidth',2);
title('(d) After 20 years', 'interpreter','latex','FontSize',12);
xlabel('assets ($a$)','interpreter','latex','FontSize',12)
ylim([ylimL ylimH]);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

legend({'HL-SSS','LL-SSS'},'Location','southeast', 'interpreter','latex','FontSize',10)



print -dpdf h55_cIRF
savefig(myfig,'h55_cIRF.fig');

print -dpdf g55_cIRF


%%
% DIRF

ylimL1 = floor(min(squeeze(IRF_hlsss_g(1:100,1,1:1/dt:40/dt)+IRF_hlsss_g(1:100,2,1:1/dt:40/dt)))*100)/100;
ylimH1 = ceil (max(squeeze(IRF_hlsss_g(1:100,1,1:1/dt:40/dt)+IRF_hlsss_g(1:100,2,1:1/dt:40/dt)))*100)/100;

ylimL2 = floor(min(squeeze(IRF_llsss_g(1:100,1,1:1/dt:40/dt)+IRF_llsss_g(1:100,2,1:1/dt:40/dt)))*100)/100;
ylimH2 = ceil (max(squeeze(IRF_llsss_g(1:100,1,1:1/dt:40/dt)+IRF_llsss_g(1:100,2,1:1/dt:40/dt)))*100)/100;

ylimL = min([ylimL1 ylimL2]);
ylimH = max([ylimH1 ylimH2]);


myfig=figure(57);
set(myfig, 'Position', [0 0 800 400])

subplot(1,2,1);
mesh([1:1/dt:40/dt]*dt,a_grid(1:100),squeeze(IRF_hlsss_g(1:100,1,1:1/dt:40/dt)+IRF_hlsss_g(1:100,2,1:1/dt:40/dt)));
title('(a) HL-SSS', 'interpreter','latex','FontSize',12);
axis([0 40 a_grid(1) a_grid(100) ylimL ylimH])
view([75, 10])
xlabel('years','interpreter','latex','FontSize',12)
ylabel('assets ($a$)','interpreter','latex','FontSize',12)

subplot(1,2,2);
mesh([1:1/dt:40/dt]*dt,a_grid(1:100),squeeze(IRF_llsss_g(1:100,1,1:1/dt:40/dt)+IRF_llsss_g(1:100,2,1:1/dt:40/dt)));
title('(b) LL-SSS', 'interpreter','latex','FontSize',12);
axis([0 40 a_grid(1) a_grid(100) ylimL ylimH])
view([75, 10])
xlabel('years','interpreter','latex','FontSize',12)
ylabel('assets ($a$)','interpreter','latex','FontSize',12)

print -dpdf h57_DIRF
savefig(myfig,'h57_DIRF.fig');

print -dpdf g57_DIRF


%% distributions, value function V, and effect of the shock on value function V

% V at llsss, after the shock

BposD =floor((irf_llsssB_B(2)-Bmin)/dB)+1;
BposU = ceil((irf_llsssB_B(2)-Bmin)/dB)+1;
wB=(B_grid(BposU)-irf_llsssB_B(2))/dB;  % weight of posD

NposD =floor((irf_llsssB_N(2)-Nmin)/dN)+1;
NposU = ceil((irf_llsssB_N(2)-Nmin)/dN)+1;
wN=(N_grid(NposU)-irf_llsssB_N(2))/dN;  % weight of posD

irf_llsss_V2 = squeeze( wB*wN*V1(:,:,BposD,NposD) + (1-wB)*wN*V1(:,:,BposU,NposD) + wB*(1-wN)*V1(:,:,BposD,NposU) + (1-wB)*(1-wN)*V1(:,:,BposU,NposU) );

% V at llsss, before the shock

BposD =floor((irf_llsssB_B(1)-Bmin)/dB)+1;
BposU = ceil((irf_llsssB_B(1)-Bmin)/dB)+1;
wB=(B_grid(BposU)-irf_llsssB_B(1))/dB;  % weight of posD

NposD =floor((irf_llsssB_N(1)-Nmin)/dN)+1;
NposU = ceil((irf_llsssB_N(1)-Nmin)/dN)+1;
wN=(N_grid(NposU)-irf_llsssB_N(1))/dN;  % weight of posD

irf_llsss_V1 = squeeze( wB*wN*V1(:,:,BposD,NposD) + (1-wB)*wN*V1(:,:,BposU,NposD) + wB*(1-wN)*V1(:,:,BposD,NposU) + (1-wB)*(1-wN)*V1(:,:,BposU,NposU) );

% V at hlsss, after the shock

BposD =floor((irf_hlsssB_B(2)-Bmin)/dB)+1;
BposU = ceil((irf_hlsssB_B(2)-Bmin)/dB)+1;
wB=(B_grid(BposU)-irf_hlsssB_B(2))/dB;  % weight of posD

NposD =floor((irf_hlsssB_N(2)-Nmin)/dN)+1;
NposU = ceil((irf_hlsssB_N(2)-Nmin)/dN)+1;
wN=(N_grid(NposU)-irf_hlsssB_N(2))/dN;  % weight of posD

irf_hlsss_V2 = squeeze( wB*wN*V1(:,:,BposD,NposD) + (1-wB)*wN*V1(:,:,BposU,NposD) + wB*(1-wN)*V1(:,:,BposD,NposU) + (1-wB)*(1-wN)*V1(:,:,BposU,NposU) );

% V at hlsss, before the shock

BposD =floor((irf_hlsssB_B(1)-Bmin)/dB)+1;
BposU = ceil((irf_hlsssB_B(1)-Bmin)/dB)+1;
wB=(B_grid(BposU)-irf_hlsssB_B(1))/dB;  % weight of posD

NposD =floor((irf_hlsssB_N(1)-Nmin)/dN)+1;
NposU = ceil((irf_hlsssB_N(1)-Nmin)/dN)+1;
wN=(N_grid(NposU)-irf_hlsssB_N(1))/dN;  % weight of posD

irf_hlsss_V1 = squeeze( wB*wN*V1(:,:,BposD,NposD) + (1-wB)*wN*V1(:,:,BposU,NposD) + wB*(1-wN)*V1(:,:,BposD,NposU) + (1-wB)*(1-wN)*V1(:,:,BposU,NposU) );






% now plot it

myfig=figure(58);
set(myfig, 'Position', [0 0 830 1000])

subplot(3,2,1);
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_llsss(2:ceil(size(g_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_hlsss(2:ceil(size(g_ss,1)/2),1)), '-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1),squeeze(g_llsss(1,1)),'o','Color',[1,0.1,0.1],'Linewidth',2)
plot(a_grid(1),squeeze(g_hlsss(1,1)),'o','Color',[0,0.5,0.1],'Linewidth',2)
title('(a) Distribution ($g$) of low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,2);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_llsss(1:ceil(size(g_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_hlsss(1:ceil(size(g_ss,1)/2),2)), '-','Color',[0,0.5,0.1],'Linewidth',2)
title('(b) Distribution ($g$) of high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'LL-SSS','HL-SSS'},'Location','northeast', 'interpreter','latex','FontSize',12)

subplot(3,2,3);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_llsss_V1(1:ceil(size(g_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_hlsss_V1(1:ceil(size(g_ss,1)/2),1)), '-','Color',[0,0.5,0.1],'Linewidth',2)
title('(c) Value function ($V$) of low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,4);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_llsss_V1(1:ceil(size(g_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_hlsss_V1(1:ceil(size(g_ss,1)/2),2)), '-','Color',[0,0.5,0.1],'Linewidth',2)
title('(d) Value function ($V$) of high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'LL-SSS','HL-SSS'},'Location','southeast', 'interpreter','latex','FontSize',12)

subplot(3,2,5);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_llsss_V2(1:ceil(size(g_ss,1)/2),1)-irf_llsss_V1(1:ceil(size(g_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_hlsss_V2(1:ceil(size(g_ss,1)/2),1)-irf_hlsss_V1(1:ceil(size(g_ss,1)/2),1)), '-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(g_ss,1)/2)),zeros(ceil(size(g_ss,1)/2),1), '-k','Linewidth',1)
title('(e) Effect of shock on $V$, low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(3,2,6);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_llsss_V2(1:ceil(size(g_ss,1)/2),2)-irf_llsss_V1(1:ceil(size(g_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(irf_hlsss_V2(1:ceil(size(g_ss,1)/2),2)-irf_hlsss_V1(1:ceil(size(g_ss,1)/2),2)), '-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(g_ss,1)/2)),zeros(ceil(size(g_ss,1)/2),1), '-k','Linewidth',1)
title('(f) Effect of shock on $V$, high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'LL-SSS','HL-SSS'},'Location','southeast', 'interpreter','latex','FontSize',12)

print -dpdf h58_g_V_DV
savefig(myfig,'h58_g_V_DV.fig');

print -dpdf g58_g_V_DV

