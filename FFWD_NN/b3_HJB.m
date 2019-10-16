% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This code solves the HJB by upwind finite differences

% VARIABLES:
% some additional 4-dimmensional matrices (the ones we want to restart at every iteration)

dVdB   = zeros(nval_a, nval_z, nval_B, nval_N);
dVdN   = zeros(nval_a, nval_z, nval_B, nval_N);
d2VddN = zeros(nval_a, nval_z, nval_B, nval_N);

dVda_f = zeros(nval_a, nval_z, nval_B, nval_N);
dVda_b = zeros(nval_a, nval_z, nval_B, nval_N);

c      = zeros(nval_a, nval_z, nval_B, nval_N);

V0     = V1;
V0_stacked = V1_stacked;

% Collections of empty sparse matrices, for building my very big A and B (here I will call them A3 and B3)
% For consistency with the paper, A1 could be called A_lm, A2 could be called A_m, A3 could be called A

for iB=1:nval_B
    for iN=1:nval_N
        A1{iB,iN}=sparse(nval_a*nval_z,nval_a*nval_z);
    end
end


A1_switch = [sparse(nval_a,nval_a),speye(nval_a)*la(1);speye(nval_a)*la(2),sparse(nval_a,nval_a)];

for iN=1:nval_N
    A2{iN}=sparse(nval_a*nval_z*nval_B,nval_a*nval_z*nval_B);
end

for iN=1:nval_N
    Xmat{iN}=sparse(nval_a*nval_z*nval_B,nval_a*nval_z*nval_B);
end

for iN=1:nval_N
    Pmat{iN}=sparse(nval_a*nval_z*nval_B,nval_a*nval_z*nval_B);
end

A3=sparse(nval_a*nval_z*nval_B*nval_N,nval_a*nval_z*nval_B*nval_N);

% HJB MAIN LOOP
for itHJB = 1:maxitHJB

    V0                  = weHJB*V1 + (1-weHJB)*V0;
    V0_stacked          = weHJB*V1_stacked + (1-weHJB)*V0_stacked;

    % forward difference
    dVda_f(1:nval_a-1,:,:,:) = (V0(2:nval_a,:,:,:)-V0(1:nval_a-1,:,:,:))/da;
    dVda_f(nval_a,:,:,:)     = (w(nval_a,:,:,:).*z(nval_a,:,:,:) + r(nval_a,:,:,:).*amax).^(-gamma);    % will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVda_b(2:nval_a,:,:,:)   = (V0(2:nval_a,:,:,:)-V0(1:nval_a-1,:,:,:))/da;
    dVda_b(1,:,:,:)          = (w(1,:,:,:).*z(1,:,:,:) + r(1,:,:,:).*amin).^(-gamma);                   % state constraint boundary condition
    
    dVdB(:,:,1:nval_B-1,:)   = (V0(:,:,2:nval_B,:)-V0(:,:,1:nval_B-1,:))/dB;                            % zero in the upper boundary
    dVdN(:,:,:,1:nval_N-1)   = (V0(:,:,:,2:nval_N)-V0(:,:,:,1:nval_N-1))/dN;                            % zero in the upper boundary
    d2VddN(:,:,:,2:nval_N-1) = (V0(:,:,:,3:nval_N)+V0(:,:,:,1:nval_N-2)-2*V0(:,:,:,2:nval_N-1))/(dN^2); % zero in the upper and lower boundaries
    
    
    I_concave = dVda_b > dVda_f;           % indicator whether value function is concave (problems arise if this is not the case) (it happens relatively often in the boundaries of a)
    
    % consumption and savings with forward difference
    cf           = (dVda_f).^(-1/gamma);
    sf           = w.*z + r.*a - cf;
    % consumption and savings with backward difference
    cb           = (dVda_b).^(-1/gamma);
    sb           = w.*z + r.*a - cb;
    % consumption and derivative of value function at non-difference
    c0           = w.*z + r.*a;
    dVda_0       = c0.^(-gamma);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If           = (sf > 0);  % positive drift --> forward difference
    Ib           = (sb < 0);  % negative drift --> backward difference
    I0           = (1-If-Ib); % at non-difference
    
    dVda_upwind  = dVda_f.*If + dVda_b.*Ib + dVda_0.*I0;
    c            = (dVda_upwind).^(-1/gamma);
    u            = (c.^(1-gamma)-1)/(1-gamma);
    
    elem_l       = zeros(nval_a, nval_z, nval_B, nval_N);
    elem_l(:,1,:,:) = la1;
    elem_l(:,2,:,:) = la2;
    
    elem_m       = alpha * Zeta * ((B+N).^alpha) - delta*(B+N) - r.*B - rhohat*N;
    elem_s       = sigma * (B+N);
    
    
    % Elements for constructing my very big sparse matrix A3
    elem_a       = (-sb.*Ib)/da;
    elem_e       = ( sf.*If)/da;
    elem_b       = -elem_a-elem_e - elem_l - PLM/dB - elem_m/dN - ((elem_s/dN).^2);
    elem_r       = ((sigma*(B+N)/dN).^2)/2;
    elem_x       = elem_m/dN + elem_r;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructing my very big sparse matrix A3 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Fill A1 collection
    
    for iB=1:nval_B
        for iN=1:nval_N
            A11          = spdiags(elem_b(:,1,iB,iN),0,nval_a,nval_a) + spdiags(elem_a(2:nval_a,1,iB,iN),-1,nval_a,nval_a) + spdiags([0;elem_e(1:nval_a-1,1,iB,iN)],1,nval_a,nval_a);
            A12          = spdiags(elem_b(:,2,iB,iN),0,nval_a,nval_a) + spdiags(elem_a(2:nval_a,2,iB,iN),-1,nval_a,nval_a) + spdiags([0;elem_e(1:nval_a-1,2,iB,iN)],1,nval_a,nval_a);
            A1{iB,iN}    = [A11,sparse(nval_a,nval_a);sparse(nval_a,nval_a),A12] + A1_switch;
        end
    end

    % Fill A2 collection
    
    jumper = nval_a*nval_z;
    for iN=1:nval_N
        for iB=1:nval_B
            A2{iN}((iB-1)*jumper+1:iB*jumper,(iB-1)*jumper+1:iB*jumper) = A1{iB,iN};
        end
    end
    for iN=1:nval_N
        for iB=1:nval_B-1
            A2{iN}((iB-1)*jumper+1:iB*jumper,iB*jumper+1:(iB+1)*jumper) = (PLM(1,1,iB,iN)/dB)*speye(nval_a*nval_z);
        end
        A2{iN}((nval_B-1)*jumper+1:nval_B*jumper,(nval_B-1)*jumper+1:nval_B*jumper)=A2{iN}((nval_B-1)*jumper+1:nval_B*jumper,(nval_B-1)*jumper+1:nval_B*jumper)+(PLM(1,1,nval_B,iN)/dB)*speye(nval_a*nval_z);
    end

    % Prepare X2 and P2 collections
    
    for iN=1:nval_N
        for iB=1:nval_B
            Xmat{iN}((iB-1)*jumper+1:iB*jumper,(iB-1)*jumper+1:iB*jumper) = elem_x(1,1,iB,iN) * speye(nval_a*nval_z);
        end
    end

    for iN=1:nval_N
        for iB=1:nval_B
            Pmat{iN}((iB-1)*jumper+1:iB*jumper,(iB-1)*jumper+1:iB*jumper) = elem_r(1,1,iB,iN) * speye(nval_a*nval_z);
        end
    end

    % Fill A3
    
    jumper = nval_a*nval_z*nval_B;
    for iN=1:nval_N
        A3((iN-1)*jumper+1:iN*jumper,(iN-1)*jumper+1:iN*jumper) = A2{iN};
    end
    for iN=1:nval_N-1
        A3(iN*jumper+1:(iN+1)*jumper,(iN-1)*jumper+1:iN*jumper) = Pmat{iN+1};
    end
    for iN=1:nval_N-1
        A3((iN-1)*jumper+1:iN*jumper,iN*jumper+1:(iN+1)*jumper) = Xmat{iN};
    end
    A3(1:jumper,1:jumper) = A3(1:jumper,1:jumper) + Pmat{1};
    A3((nval_N-1)*jumper+1:nval_N*jumper,(nval_N-1)*jumper+1:nval_N*jumper) = A3((nval_N-1)*jumper+1:nval_N*jumper,(nval_N-1)*jumper+1:nval_N*jumper) + Xmat{nval_N};
    
    % calculate B and b, find V1
    
    u_stacked = u(:);

    B3           = (rho + 1/Delta)*speye(nval_a*nval_z*nval_B*nval_N) - A3;
    
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

HJB_check = zeros(nval_a, nval_z, nval_B, nval_N);

HJB_check(:,1,:,:) = -rho*V1(:,1,:,:) + (c(:,1,:,:).^(1-gamma)-1)/(1-gamma) + (w(:,1,:,:).*z(:,1,:,:)+r(:,1,:,:).*a(:,1,:,:)-c(:,1,:,:)).*dVda_upwind(:,1,:,:) + la1*(V1(:,2,:,:)-V1(:,1,:,:)) + PLM(:,1,:,:).*dVdB(:,1,:,:) + (elem_m(:,1,:,:)).*dVdN(:,1,:,:) + (((elem_s(:,1,:,:)).^2)/2).*d2VddN(:,1,:,:);

HJB_check(:,2,:,:) = -rho*V1(:,2,:,:) + (c(:,2,:,:).^(1-gamma)-1)/(1-gamma) + (w(:,2,:,:).*z(:,2,:,:)+r(:,2,:,:).*a(:,2,:,:)-c(:,2,:,:)).*dVda_upwind(:,2,:,:) + la2*(V1(:,1,:,:)-V1(:,2,:,:)) + PLM(:,2,:,:).*dVdB(:,2,:,:) + (elem_m(:,2,:,:)).*dVdN(:,2,:,:) + (((elem_s(:,2,:,:)).^2)/2).*d2VddN(:,2,:,:);

HJB_check_stacked = zeros(nval_a*nval_z*nval_B*nval_N,1);
for iN=2:nval_N-1
    for iB=2:nval_B-1
        for iz=1:nval_z
            HJB_check_stacked(2+(iz-1)*nval_a+(iB-1)*nval_a*nval_z+(iN-1)*nval_a*nval_z*nval_B:(nval_a-1)+(iz-1)*nval_a+(iB-1)*nval_a*nval_z+(iN-1)*nval_a*nval_z*nval_B,1)=HJB_check(2:nval_a-1,iz,iB,iN);
        end
    end
end

disp('HJB_check:')
disp(max(max(abs(HJB_check_stacked))))


