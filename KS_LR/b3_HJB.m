% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This code solves the HJB by upwind finite differences

% VARIABLES:
% some additional 4-dimmensional matrices (the ones we want to restart at every iteration)

dVdK   = zeros(nval_a, nval_z, nval_K, nval_Z);
dVdZ   = zeros(nval_a, nval_z, nval_K, nval_Z);
d2VddZ = zeros(nval_a, nval_z, nval_K, nval_Z);

dVda_f = zeros(nval_a, nval_z, nval_K, nval_Z);
dVda_b = zeros(nval_a, nval_z, nval_K, nval_Z);

c      = zeros(nval_a, nval_z, nval_K, nval_Z);

u_stacked  = zeros(nval_a*nval_z*nval_K*nval_Z,1);

% Collections of empty sparse matrices, for building my very big A and B (which I will call A3 and B3)
% A1 could also be called A_lm, A2 could also be called A_m, A3 could also be called A

for iK=1:nval_K
    for iZ=1:nval_Z
        A1{iK,iZ}=sparse(nval_a*nval_z,nval_a*nval_z);
    end
end

A1_switch = [-speye(nval_a)*la(1),speye(nval_a)*la(1);speye(nval_a)*la(2),-speye(nval_a)*la(2)];

for iZ=1:nval_Z
    A2{iZ}=sparse(nval_a*nval_z*nval_K,nval_a*nval_z*nval_K);
end

A3=sparse(nval_a*nval_z*nval_K*nval_Z,nval_a*nval_z*nval_K*nval_Z);

V0=V1;

% MAIN LOOP
for itHJB = 1:maxitHJB

    V0                  = weHJB*V1 + (1-weHJB)*V0;
    V0_stacked          = V0(:);

    % forward difference
    dVda_f(1:nval_a-1,:,:,:) = (V0(2:nval_a,:,:,:)-V0(1:nval_a-1,:,:,:))/da;
    dVda_f(nval_a,:,:,:)     = (w(nval_a,:,:,:).*z(nval_a,:,:,:) + r(nval_a,:,:,:).*amax).^(-gamma);    % will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVda_b(2:nval_a,:,:,:)   = (V0(2:nval_a,:,:,:)-V0(1:nval_a-1,:,:,:))/da;
    dVda_b(1,:,:,:)          = (w(1,:,:,:).*z(1,:,:,:) + r(1,:,:,:).*amin).^(-gamma);                   % state constraint boundary condition
    
    dVdK(:,:,1:nval_K-1,:)   = (V0(:,:,2:nval_K,:)-V0(:,:,1:nval_K-1,:))/dK;                            % zero in the upper boundary
    dVdZ(:,:,:,1:nval_Z-1)   = (V0(:,:,:,2:nval_Z)-V0(:,:,:,1:nval_Z-1))/dZ;                            % zero in the upper boundary
    d2VddZ(:,:,:,2:nval_Z-1) = (V0(:,:,:,3:nval_Z)+V0(:,:,:,1:nval_Z-2)-2*V0(:,:,:,2:nval_Z-1))/(dZ^2); % zero in the upper and lower boundaries
    
    
    I_concave = dVda_b > dVda_f;           % indicator whether value function is concave (problems arise if this is not the case) (it happens relatively often in the boundaries of a)
    
    % consumption and savings with forward difference
    cf           = (dVda_f).^(-1/gamma);
    ssf          = w.*z + r.*a - cf;
    % consumption and savings with backward difference
    cb           = (dVda_b).^(-1/gamma);
    ssb          = w.*z + r.*a - cb;
    % consumption and derivative of value function at steady state
    c0           = w.*z + r.*a;
    dVda_0       = c0.^(-gamma);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If           = (ssf > 0); %positive drift --> forward difference
    Ib           = (ssb < 0); %negative drift --> backward difference
    I0           = (1-If-Ib); %at steady state
    % make sure backward difference is used at amax
    % Ib(I,:) = 1; If(I,:) = 0;
    
    dVda_upwind  = dVda_f.*If + dVda_b.*Ib + dVda_0.*I0;    %important to include third term
    c            = (dVda_upwind).^(-1/gamma);
    u            = (c.^(1-gamma)-1)/(1-gamma);
    
    % Elements for constructing my very big sparse matrix A3
    elem_a       = (-ssb.*Ib)/da;
    elem_e       = ( ssf.*If)/da;
    elem_b       = -elem_a-elem_e - PLM/dK - theta*(Zmean-Z)/dZ - ((sigma/dZ)^2);
    elem_r       = ((sigma/dZ)^2)/2;                 % this one is a scalar
    elem_x       = theta*(Zmean-Z_grid)/dZ + elem_r; % this one is a vector
    
    % Constructing my very big sparse matrix A3
    % the full collection of A1 matrices weights less than 20 MB
    % the full collection of A2 matrices weights less than 20 MB
    % A3 is a single sparse matrix, and it weights less than 25 MB

    % Fill A1 collection
    
    for iK=1:nval_K
        for iZ=1:nval_Z
            A11          = spdiags(elem_b(:,1,iK,iZ),0,nval_a,nval_a) + spdiags(elem_a(2:nval_a,1,iK,iZ),-1,nval_a,nval_a) + spdiags([0;elem_e(1:nval_a-1,1,iK,iZ)],1,nval_a,nval_a);
            A12          = spdiags(elem_b(:,2,iK,iZ),0,nval_a,nval_a) + spdiags(elem_a(2:nval_a,2,iK,iZ),-1,nval_a,nval_a) + spdiags([0;elem_e(1:nval_a-1,2,iK,iZ)],1,nval_a,nval_a);
            A1{iK,iZ}    = [A11,sparse(nval_a,nval_a);sparse(nval_a,nval_a),A12] + A1_switch;
        end
    end

    % Fill A2 collection
    
    jumper = nval_a*nval_z;
    for iZ=1:nval_Z
        for iK=1:nval_K
            A2{iZ}((iK-1)*jumper+1:iK*jumper,(iK-1)*jumper+1:iK*jumper) = A1{iK,iZ};
        end
    end
    for iZ=1:nval_Z
        for iK=1:nval_K-1
            A2{iZ}((iK-1)*jumper+1:iK*jumper,iK*jumper+1:(iK+1)*jumper) = (PLM(1,1,iK,iZ)/dK)*speye(nval_a*nval_z);
        end
        A2{iZ}((nval_K-1)*jumper+1:nval_K*jumper,(nval_K-1)*jumper+1:nval_K*jumper)=A2{iZ}((nval_K-1)*jumper+1:nval_K*jumper,(nval_K-1)*jumper+1:nval_K*jumper)+(PLM(1,1,nval_K,iZ)/dK)*speye(nval_a*nval_z);
    end

    % Fill A3
    
    jumper = nval_a*nval_z*nval_K;
    for iZ=1:nval_Z
        A3((iZ-1)*jumper+1:iZ*jumper,(iZ-1)*jumper+1:iZ*jumper) = A2{iZ};
    end
    for iZ=1:nval_Z-1
        A3(iZ*jumper+1:(iZ+1)*jumper,(iZ-1)*jumper+1:iZ*jumper) = elem_r*speye(nval_a*nval_z*nval_K);
    end
    for iZ=1:nval_Z-1
        A3((iZ-1)*jumper+1:iZ*jumper,iZ*jumper+1:(iZ+1)*jumper) = elem_x(iZ)*speye(nval_a*nval_z*nval_K);
    end
    A3((nval_Z-1)*jumper+1:nval_Z*jumper,(nval_Z-1)*jumper+1:nval_Z*jumper) = A3((nval_Z-1)*jumper+1:nval_Z*jumper,(nval_Z-1)*jumper+1:nval_Z*jumper) + elem_x(nval_Z)*speye(nval_a*nval_z*nval_K);
    A3(1:jumper,1:jumper) = A3(1:jumper,1:jumper) + elem_r*speye(nval_a*nval_z*nval_K);
    
    % calculate B and b, find V1
    
    u_stacked = u(:);

    B3           = (rho + 1/Delta)*speye(nval_a*nval_z*nval_K*nval_Z) - A3;
    
    d3           = u_stacked + V0_stacked/Delta;
    
    V1_stacked   = B3\d3; %SOLVE SYSTEM OF EQUATIONS
    
    V1=reshape(V1_stacked,size(V0));

    Vchange      = V1_stacked - V0_stacked;

    % CHECK CONVERGENCE
    dist(itHJB) = max(abs(Vchange));
    disp(dist(itHJB))
    if dist(itHJB)<critHJB
        break
    end    
end
if itHJB == maxitHJB
     disp('maximum number of iterations HJB ')
end

%%
% check that the algorithm has converged correctly

HJB_check = zeros(nval_a, nval_z, nval_K, nval_Z);

HJB_check(:,1,:,:) = -rho*V1(:,1,:,:) + (c(:,1,:,:).^(1-gamma)-1)/(1-gamma) + (w(:,1,:,:).*z(:,1,:,:)+r(:,1,:,:).*a(:,1,:,:)-c(:,1,:,:)).*dVda_upwind(:,1,:,:) + la1*(V1(:,2,:,:)-V1(:,1,:,:)) + PLM(:,1,:,:).*dVdK(:,1,:,:) + theta*(Zmean-Z(:,1,:,:)).*dVdZ(:,1,:,:) + ((sigma^2)/2)*d2VddZ(:,1,:,:);

HJB_check(:,2,:,:) = -rho*V1(:,2,:,:) + (c(:,2,:,:).^(1-gamma)-1)/(1-gamma) + (w(:,2,:,:).*z(:,2,:,:)+r(:,2,:,:).*a(:,2,:,:)-c(:,2,:,:)).*dVda_upwind(:,2,:,:) + la2*(V1(:,1,:,:)-V1(:,2,:,:)) + PLM(:,2,:,:).*dVdK(:,2,:,:) + theta*(Zmean-Z(:,2,:,:)).*dVdZ(:,2,:,:) + ((sigma^2)/2)*d2VddZ(:,2,:,:);

HJB_check_stacked = zeros(nval_a*nval_z*nval_K*nval_Z,1);
for iZ=2:nval_Z-1    % reshape could be faster, but this is clearer
    for iK=2:nval_K-1
        for iz=1:nval_z
            HJB_check_stacked(2+(iz-1)*nval_a+(iK-1)*nval_a*nval_z+(iZ-1)*nval_a*nval_z*nval_K:(nval_a-1)+(iz-1)*nval_a+(iK-1)*nval_a*nval_z+(iZ-1)*nval_a*nval_z*nval_K,1)=HJB_check(2:nval_a-1,iz,iK,iZ);
        end
    end
end

disp('HJB_check:')
disp(max(max(abs(HJB_check_stacked))))


