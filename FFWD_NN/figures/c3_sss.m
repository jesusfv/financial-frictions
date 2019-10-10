% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

sss_nval_sim=5000/dt;

tic

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first calculate some more variables for the deterministic SS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_ss = B_ss + N_ss;
r_ss = alpha * Zeta * (K_ss.^(alpha-1)) - delta - sigma2*(K_ss./N_ss);

B_ss_posD =floor((B_ss-Bmin)/dB)+1;
B_ss_posD = max(B_ss_posD,1);
B_ss_posD = min(B_ss_posD,nval_B-1);
B_ss_posU = B_ss_posD+1;
wB=(B_grid(B_ss_posU)-B_ss)/dB;

N_ss_posD =floor((N_ss-Nmin)/dN)+1;
N_ss_posD = max(N_ss_posD,1);
N_ss_posD = min(N_ss_posD,nval_N-1);
N_ss_posU = N_ss_posD+1;
wN=(N_grid(N_ss_posU)-N_ss)/dN;

c_DD=squeeze(c(:,:,B_ss_posD,N_ss_posD));
c_DU=squeeze(c(:,:,B_ss_posD,N_ss_posU));
c_UD=squeeze(c(:,:,B_ss_posU,N_ss_posD));
c_UU=squeeze(c(:,:,B_ss_posU,N_ss_posU));
c_ss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,B_ss_posD,N_ss_posD));
w_DU=squeeze(w(:,:,B_ss_posD,N_ss_posU));
w_UD=squeeze(w(:,:,B_ss_posU,N_ss_posD));
w_UU=squeeze(w(:,:,B_ss_posU,N_ss_posU));
w_ss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_ss = w_ss.*squeeze(z(:,:,1,1)) + r_ss*squeeze(a(:,:,1,1)) - c_ss;


%%
%%%%%%%%%%%%%%%%%%%%%%%
% refine HL-SSS point %
%%%%%%%%%%%%%%%%%%%%%%%

% find the point in the simulation that is closest to the HL-SSS, to start the simulation from there

mydistance = sqrt((Bsim-B_hlsss).^2 + (Nsim-N_hlsss).^2);
[mydistance_,myindex]=min(mydistance);
clear mydistance mydistance_

% now simulate without any shocks

for t=1:sss_nval_sim
    
    if t==1
        g1=squeeze(g_big(:,:,myindex));
        sss_Nsim=Nsim(myindex);
    else
        
        g0=g1;

        myA=A1{sss_BposD,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposD,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;

        sss_Nsim = sss_Nsim + ( alpha*Zeta*((sss_Bsim+sss_Nsim).^alpha) - delta*(sss_Bsim+sss_Nsim) - sss_rsim*sss_Bsim - rhohat*sss_Nsim )*dt;

    end

    sss_Bsim=sum(sum(g1.*a2*da));
    sss_Bsim=max([sss_Bsim Bmin]);
    sss_Bsim=min([sss_Bsim Bmax]);
    sss_BposD =floor((sss_Bsim-Bmin)/dB)+1;
    sss_BposU = ceil((sss_Bsim-Bmin)/dB)+1;
    wB=(B_grid(sss_BposU)-sss_Bsim)/dB;

    sss_Nsim=max([sss_Nsim Nmin]);
    sss_Nsim=min([sss_Nsim Nmax]);
    sss_NposD =floor((sss_Nsim-Nmin)/dN)+1;
    sss_NposU = ceil((sss_Nsim-Nmin)/dN)+1;
    wN=(N_grid(sss_NposU)-sss_Nsim)/dN;
    
    sss_rsim = alpha * Zeta * ((sss_Bsim+sss_Nsim).^(alpha-1)) - delta - sigma2*((sss_Bsim+sss_Nsim)./sss_Nsim);
    
end


disp('HL-SSS refinement')
disp(['B_hlsss went from ' num2str(B_hlsss) ' to ' num2str(sss_Bsim)])
disp(['N_hlsss went from ' num2str(N_hlsss) ' to ' num2str(sss_Nsim)])
disp(' ')

B_hlsss=sss_Bsim;
N_hlsss=sss_Nsim;
K_hlsss=B_hlsss+N_hlsss;
g_hlsss=g1;

r_hlsss = sss_rsim;
c_DD=squeeze(c(:,:,sss_BposD,sss_NposD));
c_DU=squeeze(c(:,:,sss_BposD,sss_NposU));
c_UD=squeeze(c(:,:,sss_BposU,sss_NposD));
c_UU=squeeze(c(:,:,sss_BposU,sss_NposU));
c_hlsss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,sss_BposD,sss_NposD));
w_DU=squeeze(w(:,:,sss_BposD,sss_NposU));
w_UD=squeeze(w(:,:,sss_BposU,sss_NposD));
w_UU=squeeze(w(:,:,sss_BposU,sss_NposU));
w_hlsss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_hlsss = w_hlsss.*squeeze(z(:,:,1,1)) + r_hlsss*squeeze(a(:,:,1,1)) - c_hlsss;


%%
%%%%%%%%%%%%%%%%%%%%%%%
% refine LL-SSS point %
%%%%%%%%%%%%%%%%%%%%%%%

% find the point in the simulation that is closest to the LL-SSS, to start the simulation from there

mydistance = sqrt((Bsim-B_llsss).^2 + (Nsim-N_llsss).^2);
[mydistance_,myindex]=min(mydistance);
clear mydistance mydistance_

% now simulate without any shocks

for t=1:sss_nval_sim
    
    if t==1
        g1=squeeze(g_big(:,:,myindex));
        sss_Nsim=Nsim(myindex);
    else
        
        g0=g1;

        myA=A1{sss_BposD,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposD,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;

        sss_Nsim = sss_Nsim + ( alpha*Zeta*((sss_Bsim+sss_Nsim).^alpha) - delta*(sss_Bsim+sss_Nsim) - sss_rsim*sss_Bsim - rhohat*sss_Nsim )*dt;

    end

    sss_Bsim=sum(sum(g1.*a2*da));
    sss_Bsim=max([sss_Bsim Bmin]);
    sss_Bsim=min([sss_Bsim Bmax]);
    sss_BposD =floor((sss_Bsim-Bmin)/dB)+1;
    sss_BposU = ceil((sss_Bsim-Bmin)/dB)+1;
    wB=(B_grid(sss_BposU)-sss_Bsim)/dB;

    sss_Nsim=max([sss_Nsim Nmin]);
    sss_Nsim=min([sss_Nsim Nmax]);
    sss_NposD =floor((sss_Nsim-Nmin)/dN)+1;
    sss_NposU = ceil((sss_Nsim-Nmin)/dN)+1;
    wN=(N_grid(sss_NposU)-sss_Nsim)/dN;
    
    sss_rsim = alpha * Zeta * ((sss_Bsim+sss_Nsim).^(alpha-1)) - delta - sigma2*((sss_Bsim+sss_Nsim)./sss_Nsim);
    
end


disp('HL-SSS refinement')
disp(['B_llsss went from ' num2str(B_llsss) ' to ' num2str(sss_Bsim)])
disp(['N_llsss went from ' num2str(N_llsss) ' to ' num2str(sss_Nsim)])
disp(' ')

B_llsss=sss_Bsim;
N_llsss=sss_Nsim;
K_llsss=B_llsss+N_llsss;
g_llsss=g1;

r_llsss = sss_rsim;
c_DD=squeeze(c(:,:,sss_BposD,sss_NposD));
c_DU=squeeze(c(:,:,sss_BposD,sss_NposU));
c_UD=squeeze(c(:,:,sss_BposU,sss_NposD));
c_UU=squeeze(c(:,:,sss_BposU,sss_NposU));
c_llsss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,sss_BposD,sss_NposD));
w_DU=squeeze(w(:,:,sss_BposD,sss_NposU));
w_UD=squeeze(w(:,:,sss_BposU,sss_NposD));
w_UU=squeeze(w(:,:,sss_BposU,sss_NposU));
w_llsss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_llsss = w_llsss.*squeeze(z(:,:,1,1)) + r_llsss*squeeze(a(:,:,1,1)) - c_llsss;


%%
%%%%%%%%%%%%%%%%%%%%
% plot SSS results %
%%%%%%%%%%%%%%%%%%%%


% distribution of assets 

myfig=figure(40);
set(myfig, 'Position', [0 0 800 400])

subplot(1,2,1);
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_llsss(2:ceil(size(g_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_hlsss(2:ceil(size(g_ss,1)/2),1)),'-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_ss   (2:ceil(size(g_ss,1)/2),1)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
plot(a_grid(1),squeeze(g_llsss(1,1)),'o','Color',[1,0.1,0.1],'Linewidth',2)
plot(a_grid(1),squeeze(g_hlsss(1,1)),'o','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1),squeeze(g_ss   (1,1)),'o','Color',[0.1,0.3,1],'Linewidth',2)
title('(a) Low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(1,2,2);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_llsss(1:ceil(size(g_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_hlsss(1:ceil(size(g_ss,1)/2),2)),'-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_ss   (1:ceil(size(g_ss,1)/2),2)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(b) High-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'LL-SSS','HL-SSS','Deterministic SS',},'Location','northeast', 'interpreter','latex','FontSize',12)

print -dpdf h40_SSS
savefig(myfig,'h40_SSS.fig');

print -dpdf g40_SSS

disp('Time refining SSS points and plotting SSS graphs:')
toc
disp(' ')

%%




disp('Gini coefficients:')

g0=ones(size(g_hlsss(:,1)))/nval_a/da; % uniform, just for testing my Gini calculations

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['test (uniform, should be 1/3):    ' num2str(g_gini) '    '  num2str(g_gini_)])

myfig = figure(41);
plot(cumsum(gini_1),gini_4)
hold on



g0 = g_hlsss(:,1)+g_hlsss(:,2);

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['HL-SSS:    ' num2str(g_gini) '    '  num2str(g_gini_)])

plot(cumsum(gini_1),gini_4)



g0 = g_llsss(:,1) + g_llsss(:,2);

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['LL-SSS:    ' num2str(g_gini) '    '  num2str(g_gini_)])

plot(cumsum(gini_1),gini_4)



g0 = g_ss(:,1) + g_ss(:,2);

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['DSS:    ' num2str(g_gini) '    '  num2str(g_gini_)])

plot(cumsum(gini_1),gini_4)
title('Lorenz curves (used for calculating Gini coefficients)', 'interpreter','latex','FontSize',14);
legend({'Uniform','HL-SSS','LL-SSS','DSS',},'Location','northwest', 'interpreter','latex','FontSize',12)


