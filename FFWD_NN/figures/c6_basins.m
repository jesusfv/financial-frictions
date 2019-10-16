% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

copyfile '../f8_KFE_sim.m';

basin_sim=800/dt;

basin_nB=15;
basin_nN=basin_nB;

N_set = Nmin:(Nmax-Nmin)/(basin_nN+1):Nmax;
N_set = N_set(2:end-1);

B_set = Bmin:(Bmax-Bmin)/(basin_nB+1):Bmax;
B_set = B_set(2:end-1);

g0_set=zeros(nval_a,nval_z,basin_nB);

g_reducer = 1.01:-0.02/(nval_a-1):0.99;
g_reducer = [g_reducer' g_reducer'];

for it_g = 1:basin_nB
    g0 = g_hlsss;
    gsum = sum(sum(g0.*a_grid.*da));
    found = 0;
    g_exp = 0;
    g_exp_=0.1;
    while found == 0
        if gsum > B_set(it_g)
            g_exp=g_exp+g_exp_;
        else
            g_exp=g_exp-g_exp_;
        end
        g0 = g_hlsss .* (g_reducer.^g_exp);
        g0 = g0 ./ sum(sum(g0*da));
        gsum = sum(sum(g0.*a_grid.*da));
        if abs(gsum-B_set(it_g))<0.001
            found = 1;
        end
    end
    g0_set(:,:,it_g)=g0;
end
    
basin_Bsim_big    = zeros(basin_sim,basin_nB*basin_nN);
basin_Nsim_big    = zeros(basin_sim,basin_nB*basin_nN);
basin_Ksim_big    = zeros(basin_sim,basin_nB*basin_nN);

KFE_inputs.each_sim = basin_sim;
KFE_inputs.nval_a   = nval_a;
KFE_inputs.nval_z   = nval_z;
KFE_inputs.A1       = A1;
KFE_inputs.a2       = a2;
KFE_inputs.dt       = dt;
KFE_inputs.da       = da;
KFE_inputs.dB       = dB;
KFE_inputs.dN       = dN;
KFE_inputs.dBB      = dBB;
KFE_inputs.dNN      = dNN;
KFE_inputs.Bmin     = Bmin;
KFE_inputs.Bmax     = Bmax;
KFE_inputs.Nmin     = Nmin;
KFE_inputs.Nmax     = Nmax;
KFE_inputs.alpha    = alpha;
KFE_inputs.Zeta     = Zeta;
KFE_inputs.delta    = delta;
KFE_inputs.sigma2   = sigma2;
KFE_inputs.rhohat   = rhohat;
KFE_inputs.B_grid   = B_grid;
KFE_inputs.N_grid   = N_grid;
KFE_inputs.sdte_shocks = zeros(basin_sim,1);

if use_parallel == 0

    disp('KFE computations will NOT use the parallel toolbox');

    simpos=1;
    for iN=1:basin_nN
        for iB=1:basin_nB

            KFE_inputs.g_ss     = g0_set(:,:,iB);
            KFE_inputs.N_ss     = N_set(iN);

            KFE_sim = f8_KFE_sim(KFE_inputs);

            basin_Bsim_big(:,simpos)=KFE_sim.Bsim;
            basin_Nsim_big(:,simpos)=KFE_sim.Nsim;

            simpos=simpos+1;
            
        end

    end
    
    clear KFE_sim KFE_inputs; 

else
    
    disp('KFE computations will use the parallel toolbox');

    myparpool = gcp();

    simpos=1;
    for iN=1:basin_nN
        for iB=1:basin_nB

            KFE_inputs.g_ss     = g0_set(:,:,iB);
            KFE_inputs.N_ss     = N_set(iN);

            f(simpos) = parfeval(myparpool,@f8_KFE_sim,1,KFE_inputs);

            simpos=simpos+1;

        end
    end

    simpos=1;
    for iN=1:basin_nN
        for iB=1:basin_nB

            KFE_sim = fetchOutputs(f(simpos));

            basin_Bsim_big(:,simpos)=KFE_sim.Bsim;
            basin_Nsim_big(:,simpos)=KFE_sim.Nsim;

            simpos=simpos+1;

        end
    end
    
    clear KFE_sim KFE_inputs;
    delete(gcp('nocreate'));
    
end

basin_B  = zeros(basin_nB*basin_nN,1);
basin_N  = zeros(basin_nB*basin_nN,1);
basin_ll = zeros(basin_nB*basin_nN,1);
for it_s = 1:basin_nB*basin_nN
    basin_B(it_s)=basin_Bsim_big(1,it_s);
    basin_N(it_s)=basin_Nsim_big(1,it_s);
    dist_hl = ((basin_Bsim_big(end,it_s)-B_hlsss).^2 + (basin_Nsim_big(end,it_s)-N_hlsss).^2).^0.5;
    dist_ll = ((basin_Bsim_big(end,it_s)-B_llsss).^2 + (basin_Nsim_big(end,it_s)-N_llsss).^2).^0.5;
    if dist_ll<dist_hl
        basin_ll(it_s)=1;
    end
end

basin_Ksim_big = basin_Bsim_big + basin_Nsim_big;

%%

myfig=figure(37);
set(myfig, 'Position', [0 0 800 800])
plot(B_hlsss,N_hlsss,'ko')
hold on
plot(B_llsss,N_llsss,'ko')
for it_s = 1:basin_nB*basin_nN
    if it_s == 122             % plot a longer line for this one, just to avoid an ugly gap in the graph
        if basin_ll(it_s)==1
            plot(basin_Bsim_big(1:1/dt:end,it_s),basin_Nsim_big(1:1/dt:end,it_s),'-','Color',[1,0.1,0.1])
        else
            plot(basin_Bsim_big(1:1/dt:end,it_s),basin_Nsim_big(1:1/dt:end,it_s),'-','Color',[0,0.5,0.1])
        end
    else                       % simulations had to be long to assure closeness to its SSS, but we plot fewer points to get a higher-quality output file
        if basin_ll(it_s)==1
            plot(basin_Bsim_big(1:1/dt:end*0.3,it_s),basin_Nsim_big(1:1/dt:end*0.3,it_s),'-','Color',[1,0.1,0.1])
        else
            plot(basin_Bsim_big(1:1/dt:end*0.3,it_s),basin_Nsim_big(1:1/dt:end*0.3,it_s),'-','Color',[0,0.5,0.1])
        end
    end
        
end
title('Convergence paths on the ($B$,$N$) space', 'interpreter','latex','FontSize',14)
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14)
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14)
text(0.85,3.0,'Basin of attraction, LL-SSS','interpreter','latex','FontSize',14, 'Color',[0.7 0 0])
text(1.00,1.5,'Basin of attraction, HL-SSS','interpreter','latex','FontSize',14, 'Color',[0 0.3 0])
xlim([Bmin Bmax]);
ylim([Nmin Nmax]);

print -dpdf h37_basins
savefig(myfig,'h37_basins.fig');

title(' ', 'interpreter','latex','FontSize',14)
print -dpdf g37_basins
