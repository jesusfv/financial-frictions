close all

myfig=figure(42);
set(myfig, 'Position', [0 0 800 800])

subplot(2,2,1);
plot(a_grid(1:ceil(size(c_ss,1)/2)),squeeze(c_llsss(1:ceil(size(c_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(c_ss,1)/2)),squeeze(c_hlsss(1:ceil(size(c_ss,1)/2),1)), '-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(c_ss,1)/2)),squeeze(c_ss   (1:ceil(size(c_ss,1)/2),1)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(a) Consumpion ($c$) of low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(2,2,2);
plot(a_grid(1:ceil(size(c_ss,1)/2)),squeeze(c_llsss(1:ceil(size(c_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(c_ss,1)/2)),squeeze(c_hlsss(1:ceil(size(c_ss,1)/2),2)), '-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(c_ss,1)/2)),squeeze(c_ss   (1:ceil(size(c_ss,1)/2),2)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(b) Consumption ($c$) of high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(2,2,3);
plot(a_grid(1:ceil(size(s_ss,1)/2)),squeeze(s_llsss(1:ceil(size(s_ss,1)/2),1)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(s_ss,1)/2)),squeeze(s_hlsss(1:ceil(size(s_ss,1)/2),1)), '-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(s_ss,1)/2)),squeeze(s_ss   (1:ceil(size(s_ss,1)/2),1)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(c) Savings ($s$) of low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(2,2,4);
plot(a_grid(1:ceil(size(s_ss,1)/2)),squeeze(s_llsss(1:ceil(size(s_ss,1)/2),2)),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(s_ss,1)/2)),squeeze(s_hlsss(1:ceil(size(s_ss,1)/2),2)),'-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(s_ss,1)/2)),squeeze(s_ss   (1:ceil(size(s_ss,1)/2),2)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(d) Savings ($s$) of high-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'LL-SSS','HL-SSS','DSS',},'Location','northeast', 'interpreter','latex','FontSize',12)

print -dpdf h42_policyfunctions_c_s
savefig(myfig,'h42_policyfunctions_c_s.fig');

print -dpdf g42_policyfunctions_c_s






