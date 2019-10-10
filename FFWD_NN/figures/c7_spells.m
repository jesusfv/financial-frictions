% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

% assign each simulated point to a basin

basin_MNR = mnrfit([basin_B basin_N],basin_ll+1);
basin_MNRfit = [ones(size(Bsim)) Bsim Nsim] * basin_MNR;

sim_zone = (basin_MNRfit >0); % one if the point belongs to the hlsss basin
clear basin_MNR basin_MNRfit

disp('percentage of periods spent on the hlsss basin instead of the llsss basin:')
mean(sim_zone)*100

figure;
plot(sim_zone);


% now see how long spells are around one sss or the other

spells = [];
spells_hlsss = [];
spells_llsss = [];
for it10=1:multi_sim
    last = (it10-1)*each_sim+delay_sim+1;
    neighborhood = sim_zone(last);
    for it11=delay_sim+2:each_sim
        now_= (it10-1)*each_sim+it11;
        if sim_zone(now_)~=neighborhood
            spells = [spells ; [neighborhood (now_-last)*dt]];
            if neighborhood == 1
                spells_hlsss=[spells_hlsss ; (now_-last)*dt];
            else
                spells_llsss=[spells_llsss ; (now_-last)*dt];
            end
            last = now_;
            neighborhood = sim_zone(last);
        end
    end
end

disp('mean duration of spells around hlsss and llsss:')
mean(spells_hlsss)
mean(spells_llsss)

%%
% plot results

myfig=figure(38);
set(myfig, 'Position', [0 0 600 400])
subplot(2,1,1)
histogram(spells_hlsss,[1 2 4 8 16 32 64 128 256 512 1024 2048 4096])
set(gca,'xscale','log') 
title(['Average duration of spells on HL-SSS basin: ' num2str(mean(spells_hlsss)) ' years'], 'interpreter','latex','FontSize',12);
ylabel('number of spells', 'interpreter','latex','FontSize',12);
xlabel('years', 'interpreter','latex','FontSize',12);
subplot(2,1,2)
histogram(spells_llsss,[1 2 4 8 16 32 64 128 256 512 1024 2048 4096])
set(gca,'xscale','log') 
title(['Average duration of spells on LL-SSS basin: ' num2str(mean(spells_llsss)) ' years'], 'interpreter','latex','FontSize',12);
ylabel('number of spells', 'interpreter','latex','FontSize',12);
xlabel('years', 'interpreter','latex','FontSize',12);

print -dpdf h38_spells
savefig(myfig,'h38_spells.fig');

print -dpdf g38_spells

%%

wsim = (1-alpha) * Zeta * (Bsim+Nsim).^alpha;

Bsim_hlsss = Bsim .* (1./sim_zone);
Nsim_hlsss = Nsim .* (1./sim_zone);
Bsim_llsss= Bsim .* (1./(1-sim_zone));
Nsim_llsss= Nsim .* (1./(1-sim_zone));


myfig=figure(36);
set(myfig, 'Position', [0 0 800 800])
plot(Bsim_hlsss(1:used_sim),Nsim_hlsss(1:used_sim),'-','Color',[0,0.5,0.1],'linewidth',1);
hold on;
plot(Bsim_llsss(1:used_sim),Nsim_llsss(1:used_sim),'-','Color',[1,0.1,0.1],'linewidth',1);
title('Simulation', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
text(1.0,3.0,'Basin of attraction, LL-SSS','interpreter','latex','FontSize',14, 'Color',[0.7 0 0])
text(1.7,2.5,'Basin of attraction, HL-SSS','interpreter','latex','FontSize',14, 'Color',[0 0.3 0])
xlim([Bmin Bmax])
ylim([Nmin Nmax])

print -dpdf h36_sim
savefig(myfig,'h36_sim.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g36_sim

%%

Ysim_hlsss = Ysim(sim_zone);
rsim_hlsss = rsim(sim_zone);
wsim_hlsss = wsim(sim_zone);
sim_zone_=not(sim_zone);
Ysim_llsss= Ysim(sim_zone_);
rsim_llsss= rsim(sim_zone_);
wsim_llsss= wsim(sim_zone_);

myfig=figure(39);
set(myfig, 'Position', [0 0 800 800])

subplot(2,2,1:2)
[y_ax,x_ax] = hist(Ysim_hlsss,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'-','Color',[0,0.5,0.1],'Linewidth',2)
hold on
[y_ax,x_ax] = hist(Ysim_llsss,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'--','Color',[1,0.1,0.1],'Linewidth',2)
grid
ylabel('frequency', 'interpreter','latex','FontSize',14);
title('(a) Output ($Y$)', 'interpreter','latex','FontSize',14);
legend({'HL-SSS basin','LL-SSS basin'},'Location','northwest', 'interpreter','latex','FontSize',12)

subplot(2,2,3)
[y_ax,x_ax] = hist(rsim_hlsss,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'-','Color',[0,0.5,0.1],'Linewidth',2)
hold on
[y_ax,x_ax] = hist(rsim_llsss,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'--','Color',[1,0.1,0.1],'Linewidth',2)
grid
ylabel('frequency', 'interpreter','latex','FontSize',14);
title('(b) Interest rate ($r$)', 'interpreter','latex','FontSize',14);

subplot(2,2,4)
[y_ax,x_ax] = hist(wsim_hlsss,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'-','Color',[0,0.5,0.1],'Linewidth',2)
hold on
[y_ax,x_ax] = hist(wsim_llsss,60);
dx = x_ax(2)-x_ax(1);
y_ax = y_ax / sum(y_ax*dx);
plot(x_ax,y_ax,'--','Color',[1,0.1,0.1],'Linewidth',2)
grid
ylabel('frequency', 'interpreter','latex','FontSize',14);
title('(c) Wages ($w$)', 'interpreter','latex','FontSize',14);

print -dpdf h41_sim_distr
savefig(myfig,'h41_sim_distr.fig');

print -dpdf g41_sim_distr

%%

disp(' ')
disp(['mean of Y:           ' num2str(mean(Ysim))])
disp(['mean of Y at HL-SSS: ' num2str(mean(Ysim_hlsss))])
disp(['mean of Y at LL-SSS: ' num2str(mean(Ysim_llsss))])
disp(' ')
disp(['mean of r:           ' num2str(mean(rsim))])
disp(['mean of r at HL-SSS: ' num2str(mean(rsim_hlsss))])
disp(['mean of r at LL-SSS: ' num2str(mean(rsim_llsss))])
disp(' ')
disp(['mean of w:           ' num2str(mean(wsim))])
disp(['mean of w at HL-SSS: ' num2str(mean(wsim_hlsss))])
disp(['mean of w at LL-SSS: ' num2str(mean(wsim_llsss))])
disp(' ')

disp(' ')
disp(['std of Y:           ' num2str(std(Ysim))])
disp(['std of Y at HL-SSS: ' num2str(std(Ysim_hlsss))])
disp(['std of Y at LL-SSS: ' num2str(std(Ysim_llsss))])
disp(' ')
disp(['std of r:           ' num2str(std(rsim))])
disp(['std of r at HL-SSS: ' num2str(std(rsim_hlsss))])
disp(['std of r at LL-SSS: ' num2str(std(rsim_llsss))])
disp(' ')
disp(['std of w:           ' num2str(std(wsim))])
disp(['std of w at HL-SSS: ' num2str(std(wsim_hlsss))])
disp(['std of w at LL-SSS: ' num2str(std(wsim_llsss))])
disp(' ')

disp(' ')
disp(['skewness of Y:           ' num2str(skewness(Ysim))])
disp(['skewness of Y at HL-SSS: ' num2str(skewness(Ysim_hlsss))])
disp(['skewness of Y at LL-SSS: ' num2str(skewness(Ysim_llsss))])
disp(' ')
disp(['skewness of r:           ' num2str(skewness(rsim))])
disp(['skewness of r at HL-SSS: ' num2str(skewness(rsim_hlsss))])
disp(['skewness of r at LL-SSS: ' num2str(skewness(rsim_llsss))])
disp(' ')
disp(['skewness of w:           ' num2str(skewness(wsim))])
disp(['skewness of w at HL-SSS: ' num2str(skewness(wsim_hlsss))])
disp(['skewness of w at LL-SSS: ' num2str(skewness(wsim_llsss))])
disp(' ')

disp(' ')
disp(['kurtosis of Y:           ' num2str(kurtosis(Ysim))])
disp(['kurtosis of Y at HL-SSS: ' num2str(kurtosis(Ysim_hlsss))])
disp(['kurtosis of Y at LL-SSS: ' num2str(kurtosis(Ysim_llsss))])
disp(' ')
disp(['kurtosis of r:           ' num2str(kurtosis(rsim))])
disp(['kurtosis of r at HL-SSS: ' num2str(kurtosis(rsim_hlsss))])
disp(['kurtosis of r at LL-SSS: ' num2str(kurtosis(rsim_llsss))])
disp(' ')
disp(['kurtosis of w:           ' num2str(kurtosis(wsim))])
disp(['kurtosis of w at HL-SSS: ' num2str(kurtosis(wsim_hlsss))])
disp(['kurtosis of w at LL-SSS: ' num2str(kurtosis(wsim_llsss))])
disp(' ')

clear wsim Bsim_hlsss Nsim_hlsss Ksim_hlsss Ysim_hlsss rsim_hlsss wsim_hlsss Bsim_llsss Nsim_llsss Ksim_llsss Ysim_llsss rsim_llsss wsim_llsss sim_zone_;
