% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

sss_nval_sim=50000/dt;

sss_Ksim   = zeros(sss_nval_sim,1);

tic

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate from deterministic SS with no shocks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that the center position in Z_grid corresponds to Zmean

Zpos=(nval_Z+1)/2;
if Z_grid(Zpos)~=Zmean
    aslahskljhskjahaksjhfda % Just so it stops if it ever arrives here
end

% now simulate with no shocks

for t=1:sss_nval_sim
    
    if t==1
        g1=g_ss;
    else

        g0=g1;
        
        % Interpolation step to compute next period's distribution
        % The following block computes the transition g(t+1) = (I - A*dt)^(-1)*g(t) for the two closest grid points

        myA        = A1{sss_KposD,Zpos};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1D        = reshape(g1_stacked,size(g0));

        myA        = A1{sss_KposU,Zpos};
        myA        = (speye(nval_a*nval_z)-myA'*dt);
        g1_stacked = myA\g0(:);
        g1_stacked = g1_stacked/sum(g1_stacked*da);
        g1U        = reshape(g1_stacked,size(g0));

        % The we compute distribution next period, averaging those two results
        g1 = sss_wK*g1D + (1-sss_wK)*g1U;

    end
    
    sss_Ksim(t)    =sum(sum(g1.*a2*da));
    
    sss_Ksim(t)  =   max([sss_Ksim(t) Kmin+0.000001]);
    sss_Ksim(t)  =   min([sss_Ksim(t) Kmax-0.000001]);
    sss_KposD    = floor((sss_Ksim(t)-Kmin)/dK)+1;
    sss_KposU    =  ceil((sss_Ksim(t)-Kmin)/dK)+1;
    sss_wK=(K_grid(sss_KposU)-sss_Ksim(t))/dK;         % weight of KposD

end

disp('SSS refinement')
disp(['K_sss went from ' num2str(K_sss) ' to ' num2str(sss_Ksim(end))])
disp(' ')

K_sss=sss_Ksim(end);
g_sss=g1;
Z_sss=Zmean;

V_D=squeeze(V1(:,:,sss_KposD,Zpos));
V_U=squeeze(V1(:,:,sss_KposU,Zpos));
V1_sss = wK*V_D + (1-wK)*V_U;


%%
%%%%%%%%%%%%%%%%%%%%
% plot SSS results %
%%%%%%%%%%%%%%%%%%%%

% distribution of assets 

myfig=figure(40);
set(myfig, 'Position', [0 0 800 400])

subplot(1,2,1);
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_sss(2:ceil(size(g_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_ss (2:ceil(size(g_ss,1)/2),1)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
plot(a_grid(1),squeeze(g_sss(1,1)),'o','Color',[1,0.1,0.1],'Linewidth',2)
plot(a_grid(1),squeeze(g_ss (1,1)),'o','Color',[0.1,0.3,1],'Linewidth',2)
title('(a) Low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(1,2,2);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_sss(1:ceil(size(g_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_ss (1:ceil(size(g_ss,1)/2),2)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(b) High-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'SSS','Deterministic SS',},'Location','northeast', 'interpreter','latex','FontSize',12)

print -dpdf h40_SSS
savefig(myfig,'h40_SSS.fig');

print -dpdf g40_SSS


% evolution of K from SS to SSS

myfig=figure(41);
set(myfig, 'Position', [0 0 600 400])

plot(sss_Ksim)
title('Evolution of $K$ from DSS to SSS', 'interpreter','latex','FontSize',14);
xlabel('time (months)', 'interpreter','latex','FontSize',14);
grid

print -dpdf h41_K_SSS
savefig(myfig,'h41_K_SSS.fig');


disp('Time refining SSS points and plotting SSS graphs:')
toc
disp(' ')

