% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution


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











Bsim_used=[];
Nsim_used=[];
nsim_used=[];
for it_sim=1:multi_sim
    Bsim_used =[Bsim_used ; Bsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)];
    Nsim_used =[Nsim_used ; Nsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)];
    nsim_used =[nsim_used  (it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim];
end
nsim_used=nsim_used';
lvrg_sim= (Bsim_used+Nsim_used)./Nsim_used;

[lvrg_sort,lvrg_index]=sort(lvrg_sim);

sim05_=lvrg_index(round(0.05*size(lvrg_index,1)));
sim05 =nsim_used(sim05_);
g05 = squeeze(g_big(:,:,sim05));
B05 = Bsim(sim05);
N05 = Nsim(sim05);

for it1=2:19
    eval(['sim' num2str(5*it1) '_=lvrg_index(round(' num2str(0.05*it1) '*size(lvrg_index,1)));']);
    eval(['sim' num2str(5*it1) ' =nsim_used(sim' num2str(5*it1) '_);']);
    eval(['g'   num2str(5*it1) ' =squeeze(g_big(:,:,sim' num2str(5*it1) '));']);
    eval(['B'   num2str(5*it1) ' =Bsim(sim' num2str(5*it1) ');']);
    eval(['N'   num2str(5*it1) ' =Nsim(sim' num2str(5*it1) ');']);
end











sss_nval_sim=1000/dt;

sss_Bsim    = zeros(sss_nval_sim,1);
sss_BposD   = zeros(sss_nval_sim,1);
sss_BposU   = zeros(sss_nval_sim,1);

sss_Nsim    = zeros(sss_nval_sim,1);
sss_NposU   = zeros(sss_nval_sim,1);
sss_NposD   = zeros(sss_nval_sim,1);

sss_rsim    = zeros(sss_nval_sim,1);

sss_g_big   = zeros(nval_a,nval_z,sss_nval_sim+1);
sss_g_dif   = zeros(sss_nval_sim,1);


for t=1:sss_nval_sim
    
    if t==1
        g1=g_ss;
        sss_Nsim(1,1)=N_ss;
        sss_g_big(:,:,t)=g_ss;
    else
        
        g0=g1;

        myA=A1{sss_BposD(t-1,1),sss_NposD(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposD(t-1,1),sss_NposU(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU(t-1,1),sss_NposD(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU(t-1,1),sss_NposU(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;

        sss_Nsim(t) = sss_Nsim(t-1,1) + ( alpha*Zeta*((sss_Bsim(t-1,1)+sss_Nsim(t-1,1)).^alpha) - delta*(sss_Bsim(t-1,1)+sss_Nsim(t-1,1)) - sss_rsim(t-1,1)*sss_Bsim(t-1,1) - rhohat*sss_Nsim(t-1,1) )*dt;

    end

    sss_g_big(:,:,t+1)=g1;
    
    sss_Bsim(t)=sum(sum(g1.*a2*da));
    sss_Bsim(t)=max([sss_Bsim(t) Bmin]);
    sss_Bsim(t)=min([sss_Bsim(t) Bmax]);
    sss_BposD(t) =floor((sss_Bsim(t)-Bmin)/dB)+1;
    sss_BposU(t) = ceil((sss_Bsim(t)-Bmin)/dB)+1;
    wB=(B_grid(sss_BposU(t))-sss_Bsim(t))/dB;

    sss_Nsim(t)=max([sss_Nsim(t) Nmin]);
    sss_Nsim(t)=min([sss_Nsim(t) Nmax]);
    sss_NposD(t) =floor((sss_Nsim(t)-Nmin)/dN)+1;
    sss_NposU(t) = ceil((sss_Nsim(t)-Nmin)/dN)+1;
    wN=(N_grid(sss_NposU(t))-sss_Nsim(t))/dN;
    
    BBposD(t) =floor((sss_Bsim(t)-Bmin)/dBB)+1;
    BBposU(t) = ceil((sss_Bsim(t)-Bmin)/dBB)+1;
    NNposD(t) =floor((sss_Nsim(t)-Nmin)/dNN)+1;
    NNposU(t) = ceil((sss_Nsim(t)-Nmin)/dNN)+1;

    sss_rsim(t,1) = alpha * Zeta * ((sss_Bsim(t,1)+sss_Nsim(t,1)).^(alpha-1)) - delta - sigma2*((sss_Bsim(t,1)+sss_Nsim(t,1))./sss_Nsim(t,1));
    
    if t==1
        r_ss = sss_rsim(t,1);
        c_DD=squeeze(c(:,:,sss_BposD(t,1),sss_NposD(t)));
        c_DU=squeeze(c(:,:,sss_BposD(t,1),sss_NposU(t)));
        c_UD=squeeze(c(:,:,sss_BposU(t,1),sss_NposD(t)));
        c_UU=squeeze(c(:,:,sss_BposU(t,1),sss_NposU(t)));
        c_ss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
        w_DD=squeeze(w(:,:,sss_BposD(t,1),sss_NposD(t)));
        w_DU=squeeze(w(:,:,sss_BposD(t,1),sss_NposU(t)));
        w_UD=squeeze(w(:,:,sss_BposU(t,1),sss_NposD(t)));
        w_UU=squeeze(w(:,:,sss_BposU(t,1),sss_NposU(t)));
        w_ss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
        s_ss = w_ss.*squeeze(z(:,:,1,1)) + r_ss*squeeze(a(:,:,1,1)) - c_ss;
    end

end


sss_Ksim=sss_Bsim+sss_Nsim;

sss_g_big_average = squeeze(sss_g_big(:,1,:) + sss_g_big(:,2,:));

for t=2:sss_nval_sim
    %sss_g_dif(t,1)=(sum(sum((sss_g_big_average(t+1)-sss_g_big_average(t)).^2))).^0.5;
    sss_g_dif(t,1) =(sum(sum((sss_g_big(:,:,t)-sss_g_big(:,:,t-1)).^2))).^0.5;
end
    
B_sss=sss_Bsim(sss_nval_sim);
N_sss=sss_Nsim(sss_nval_sim);
K_sss=B_sss+N_sss;
g_sss=squeeze(sss_g_big(:,:,sss_nval_sim));

r_sss = sss_rsim(t-1,1);
c_DD=squeeze(c(:,:,sss_BposD(t-1,1),sss_NposD(t-1)));
c_DU=squeeze(c(:,:,sss_BposD(t-1,1),sss_NposU(t-1)));
c_UD=squeeze(c(:,:,sss_BposU(t-1,1),sss_NposD(t-1)));
c_UU=squeeze(c(:,:,sss_BposU(t-1,1),sss_NposU(t-1)));
c_sss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,sss_BposD(t-1,1),sss_NposD(t-1)));
w_DU=squeeze(w(:,:,sss_BposD(t-1,1),sss_NposU(t-1)));
w_UD=squeeze(w(:,:,sss_BposU(t-1,1),sss_NposD(t-1)));
w_UU=squeeze(w(:,:,sss_BposU(t-1,1),sss_NposU(t-1)));
w_sss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_sss = w_sss.*squeeze(z(:,:,1,1)) + r_sss*squeeze(a(:,:,1,1)) - c_sss;















mywidth=0.01;

B75_Nmin=inf;
B75_Nmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Bsim(t)-B75)<mywidth
        B75_Nmin=min([B75_Nmin Nsim(t)]);
        B75_Nmax=max([B75_Nmax Nsim(t)]);
    end
end

B_sss_Nmin=inf;
B_sss_Nmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Bsim(t)-B_sss)<mywidth
        B_sss_Nmin=min([B_sss_Nmin Nsim(t)]);
        B_sss_Nmax=max([B_sss_Nmax Nsim(t)]);
    end
end

B25_Nmin=inf;
B25_Nmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Bsim(t)-B25)<mywidth
        B25_Nmin=min([B25_Nmin Nsim(t)]);
        B25_Nmax=max([B25_Nmax Nsim(t)]);
    end
end

N75_Bmin=inf;
N75_Bmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Nsim(t)-N75)<mywidth
        N75_Bmin=min([N75_Bmin Bsim(t)]);
        N75_Bmax=max([N75_Bmax Bsim(t)]);
    end
end

N_sss_Bmin=inf;
N_sss_Bmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Nsim(t)-N_sss)<mywidth
        N_sss_Bmin=min([N_sss_Bmin Bsim(t)]);
        N_sss_Bmax=max([N_sss_Bmax Bsim(t)]);
    end
end

N25_Bmin=inf;
N25_Bmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Nsim(t)-N25)<mywidth
        N25_Bmin=min([N25_Bmin Bsim(t)]);
        N25_Bmax=max([N25_Bmax Bsim(t)]);
    end
end

