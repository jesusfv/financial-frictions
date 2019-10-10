function results = b3_HJB_stationary(parameters,K)
% This codes computes the HJB equation in a Bewley economy with nominal
% long-term debt
% Developped by Nuño (2016) based on codes by Ben Moll 

%PARAMETERS
alpha         = parameters.alpha;
delta         = parameters.delta;
gamma         = parameters.gamma;            % CRRA utility with parameter s
rho           = parameters.rho;              % discount rate
z1            = parameters.z1;               % labor productivity
z2            = parameters.z2;
la1           = parameters.la1;              % transition prob
la2           = parameters.la2;
I             = parameters.I;                % numer of points
amin          = parameters.amin;             % borrowing constraint
amax          = parameters.amax;             % max value of assets
maxit         = parameters.maxit;            % max number of iterations
crit          = parameters.crit;             % convergence ctrit 
Delta         = parameters.Delta;            % Delta in HJB

%VARIABLES
a   = linspace(amin,amax,I)';               % assets 1D  - I elements
da  = (amax-amin)/(I-1);     
z   = [z1,z2];                              % product 1D - 2 elements
la  = [la1,la2];    
aa  = [a,a];                                % assets 2D
zz  = ones(I,1)*z;                          % product 2D
dVf = zeros(I,2);
dVb = zeros(I,2);
c   = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)]; 

% Interest rates and wages
r = alpha * K^(alpha-1) - delta;  
w = (1-alpha) * K^alpha;


%INITIAL GUESS
v0(:,1) = (w.*z(1) + r.*a).^(1-gamma)/(1-gamma)/rho;
v0(:,2) = (w.*z(2) + r.*a).^(1-gamma)/(1-gamma)/rho;

v = v0;

%MAIN LOOP
for n = 1:maxit
    V            = v;
    V_n(:,:,n)   = V;                              % keeps track of the value function
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:)     = (w.*z + r.*amax).^(-gamma);    % will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:)     = (w.*z + r.*amin).^(-gamma);    % state constraint boundary condition
    
    I_concave = dVb > dVf;                         % indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf           = (dVf).^(-1/gamma);
    ssf          = w.*zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb           = (dVb).^(-1/gamma);
    ssb          = w.*zz + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0           = w.*zz + r.*aa;
    dV0          = c0.^(-gamma);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If           = ssf > 0; %positive drift --> forward difference
    Ib           = ssb < 0; %negative drift --> backward difference
    I0           = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    
    
    dV_Upwind    = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c            = (dV_Upwind).^(-1/gamma);
    u            = (c.^(1-gamma)-1)/(1-gamma);

    
    %CONSTRUCT MATRIX A
    elem_a            = ( -ssb.*Ib          )/da;  % - min(ssb,0)/da;
    elem_e            = (            ssf.*If)/da;  %   max(ssf,0)/da;
    elem_b            = -elem_a-elem_e;
%   elem_b            = ( -ssf.*If + ssb.*Ib)/da;  % - max(ssf,0)/da + min(ssb,0)/da;
    
    A1           = spdiags(elem_b(:,1),0,I,I) + spdiags(elem_a(2:I,1),-1,I,I) + spdiags([0;elem_e(1:I-1,1)],1,I,I);
    A2           = spdiags(elem_b(:,2),0,I,I) + spdiags(elem_a(2:I,2),-1,I,I) + spdiags([0;elem_e(1:I-1,2)],1,I,I);
    A            = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B            = (rho + 1/Delta)*speye(2*I) - A;
    
    u_stacked    = [u(:,1);u(:,2)];
    V_stacked    = [V(:,1);V(:,2)];
    
    b            = u_stacked + V_stacked/Delta;
    V_stacked    = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V            = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange      = V - v;
    v            = V;

    % CHECK CONVERGENCE
    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        %disp('Value Function Converged, Iteration = ')
        %disp(n)
        break
    end    
end
if n == maxit
     disp('maximum number of iterations HJB ')
end

% RESULTS
results.a  = a;
results.aa = aa;
results.z  = z;
results.zz = zz;
results.V  = V;
results.A  = A;
results.c  = c;
results.r  = r;
results.w  = w;
results.da = da;
results.dV = dV_Upwind;
