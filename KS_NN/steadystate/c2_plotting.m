%PLOTTING
% Color2  = [0 0 1]; %blue
% Color1  = [0 0.4 0.1]; %green
Color2  = [1 0 0]; %blue
Color1  = [0 0 0]; %green
style1 = '-';
style2 = '--';
figure(2)
subplot(221)
a = results.a;
plot(a/Y,results.V(:,1),'Color',Color1,'LineStyle',style1,'LineWidth',2)
hold on
plot(a/Y,results.V(:,2),'Color',Color2,'LineStyle',style2,'LineWidth',2)

alimmax = amax/Y*0.8;

%grid
xlabel('Assets, $a$','interpreter','latex','FontSize',16)
title('(a) Value function, $v(a)$','interpreter','latex','FontSize',16)
xlim([amin/Y alimmax])

subplot(222)
plot(a/Y,results.c(:,1),'Color',Color1,'LineStyle',style1,'LineWidth',2)
hold on
plot(a/Y,results.c(:,2),'Color',Color2,'LineStyle',style2,'LineWidth',2)
%grid
xlabel('Assets, $a$','interpreter','latex','FontSize',16)
title('(b) Consumption, $c(a)$','interpreter','latex','FontSize',16)
xlim([amin/Y alimmax])

subplot(223)
s  =  results.aa*results.r + results.w*results.zz-results.c;
plot(a/Y,s(:,1),'Color',Color1,'LineStyle',style1,'LineWidth',2)
hold on
plot(a/Y,s(:,2),'Color',Color2,'LineStyle',style2,'LineWidth',2)
%grid
xlabel('Assets, $a$','interpreter','latex','FontSize',16)
title('(c) Drift, $s(a)$','interpreter','latex','FontSize',16)
xlim([amin/Y alimmax])

subplot(224)
plot(a/Y,distribution(:,1),'Color',Color1,'LineStyle',style1,'LineWidth',2)
hold on
plot(a/Y,distribution(:,2),'Color',Color2,'LineStyle',style2,'LineWidth',2)
%grid
xlabel('Assets, $a$','interpreter','latex','FontSize',16)
title('(d) Distribution, $f(a)$','interpreter','latex','FontSize',16)
xlim([amin/Y alimmax])
legend('State 1','State 2')






