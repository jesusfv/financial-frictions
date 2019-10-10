function distribution = b4_KFE_stationary(parameters,results)
% This codes computes the KFE equation in a Bewley economy with nominal
% long-term debt
% Developped by Nuño and Thomas (2015) based on codes by Ben Moll 
I    = parameters.I;                % numer of points
da   = results.da;         
AT   = (results.A)';
b    = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix     = 4;
b(i_fix)  =.1;
row       = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg          = AT\b;
g_sum       = gg'*ones(2*I,1)*da;
gg          = gg./g_sum;

g           = [gg(1:I),gg(I+1:2*I)];

distribution = g;
