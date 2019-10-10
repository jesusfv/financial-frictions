% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

%MOMENTS
range = find(results.a<0,1,'last');
S = sum(sum(results.aa.*distribution *results.da));
C= sum(sum(results.c.*distribution *results.da));
L = sum(sum(results.zz.*distribution *results.da));
Y = S^alpha;
disp('Total capital (% GDP)')
disp (S/Y)
disp('Total labor')
disp (L)
disp('Total output')
disp (Y)
disp('Interest rate')
disp(results.r)

