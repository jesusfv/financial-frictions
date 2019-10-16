% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This script estimates the PLM based on simulated data

clear NN_all loss_final;

% prepare our observations (using only the observations after the initial pre-heat period of each simulation)

Y=[];
X1=[];
X2=[];

% Loop over several simulations

for it_sim=1:multi_sim
    Y = [Y  ; (Bsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)-Bsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1))/dt]; % dB : Growth rate of B
    X1= [X1 ; Bsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1)]; % B
    X2= [X2 ; Nsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1)]; % N
end

X = [X1 X2];

% PRE-PROCESSING: grid approximation, with linear regresion on each grid knot. This reduces the dimensionality of the problem avoiding unnnecesary noise

PLM_IP_finegrid = zeros(nval_BB,nval_NN)*NaN;

for iB=1:nval_BB
    for iN=1:nval_NN
        Y_=[];
        X1_=[];
        X2_=[];
        for it_sim=1:size(Y,1)
            if abs(X1(it_sim)-BB_grid(iB))<dBB/2
                if abs(X2(it_sim)-NN_grid(iN))<dNN/2          % If the simulation falls in a particular grid point
                    Y_=[Y_;Y(it_sim)];
                    X1_=[X1_;X1(it_sim)];
                    X2_=[X2_;X2(it_sim)];
                end
            end
        end
        if size(X1_,1)>5
            X0_=ones(size(X1_));
            X = [X0_ X1_ X2_];
            beta = (X'*X)^-1*X'*Y_;
            PLM_IP_finegrid(iB,iN) = [1 BB_grid(iB) NN_grid(iN)]*beta;  % We evaluate only the grid point
        end
    end
end

% We will use the knots in that grid of linear regression approximations to train the neural network
% (this helps because feeding less noisy data makes the NN training work faster, and takes random noise away from the full algorithm)

Y_=[];
X1_=[];
X2_=[];
for iB=1:nval_BB
    for iN=1:nval_NN
        if not(isnan(PLM_IP_finegrid(iB,iN)))
            Y_=[Y_;PLM_IP_finegrid(iB,iN)];
            X1_=[X1_;BB_grid_2D(iB,iN)];
            X2_=[X2_;NN_grid_2D(iB,iN)];
        end
    end
end

% ALTERNATIVE MACHINE LEARNING ALGORITHM (ROBUSTNESS):  an alternative to the neural network would be this built-in interpolant (IP)
% it's very fast but linear extrapolation on two dimmensions generates ugly ridges
F = scatteredInterpolant(X1_,X2_,Y_,'natural','linear');
Y_fit_IP = F(X1,X2);

% ----------------------------------
% ESTIMATION USING THE NEURAL NETWORK

% STEP 1: input normalization - sorry I wrote this in such a convoluted way, it's simpler than it looks like

Yr1=min(Y);
Yr2=max(Y);
Yr3=(Yr1+Yr2)/2; % Mean point (max  + min)/2
Yr4=Yr2-Yr3;     % Half of the interval: Yr4 = Yr2 - Yr3 = (max-min) /2 
X1r1=min(X1);
X1r2=max(X1);
X1r3=(X1r1+X1r2)/2;
X1r4=X1r2-X1r3;
X2r1=min(X2);
X2r2=max(X2);
X2r3=(X2r1+X2r2)/2;
X2r4=X2r2-X2r3;

Yr4  = Yr4 /reglimY;  % We rescale the width
X1r4 = X1r4/reglimX;
X2r4 = X2r4/reglimX;

Yr =(Y_ -Yr3 )/Yr4;   % (Y-mean)/interval_width_rescalated
X1r=(X1_-X1r3)/X1r4;
X2r=(X2_-X2r3)/X2r4;

Xr = [X1r X2r];

progress            = zeros(NN_iters,NN_starts);
learning_speed_full = zeros(NN_iters,NN_starts);
mynoise_            = mynoise;

% STEP 2: Estimation of the neural network coefficients

it1_=0;
for it0=1:NN_starts   % We run NN_starts networks
    learning_speed_=learning_speed;
    if NN_starts>1
        NN0=struct( 'b1',zeros(1,network_width), 'b2', zeros(1,1), 'w1', randn(size(Xr,2),network_width)*mynoise_, 'w2', randn(network_width,1)*mynoise_ ); % Random initial parameteres, only for the first iteration of the algorithm
        bw1=[NN0.b1 ; NN0.w1];
        NN0_flat = zeros((size(Xr,2)+2)*network_width+1,1);                % NN parameters in tensor form
        NN0_flat(1:(1+size(Xr,2))*network_width)=bw1(:);
        NN0_flat(end-network_width:end-1)=NN0.w2(:);
        NN0_flat(end)=NN0.b2;
        clear NN0 bw1;
    else
        NN0_flat=NN1_flat;   % if NN_starts=1 (which happens on all iterations except the first one) we initialize weights at the values from the previous NN
    end

    NN_flat=NN0_flat;
    
    for it1=1:NN_iters
        % use full batch instead of minibatch (possible because we're training with the grid knots; takes random noise away from the full algorithm)
        Yrows=Yr;
        Xrows=Xr;
        
        b7_PLM_iter  % This script improves the coefficient each iteration
        
        NN_flat=NN_flat_1;
        progress(it1,it0)=(fval_1^0.5)*Yr4;       % We store the loss function (de-normalized) at each step of the procedure
        learning_speed_full(it1,it0)=jumpsize1;   % We store the learning rate we used
    end
    
    NN_all(:,it0)=NN_flat;
    loss_final(it0,1)=f1_NN_loss(Yr,Xr,network_width,NN_flat,0); % Calculate loss function for the whole sample (different only if we used a minibatch) and store it (ignore regularization term by settimg lambda=0 just here)
    loss_final(it0,1)=(loss_final(it0,1).^0.5)*Yr4;              % We store it in form of de-normalized RMSE

end

[Y_RMSE_NN1,minloc]=min(loss_final); % We pick the NN that provided the minimum loss over the evaluation sample

NN1_flat = NN_all(:,minloc);

Yr=(Y-Yr3)/Yr4;      % (Y-mean)/interval_width_rescalated % now using all the simulated data
X1r=(X1-X1r3)/X1r4;
X2r=(X2-X2r3)/X2r4;

Xr = [X1r X2r];

% Given the optimal network, we simulate the outputs to compute the error
Yr_fit = f2_NN_eval(Xr,network_width,NN1_flat);
Y_fit = (Yr_fit*Yr4)+Yr3;
Y_RMSE_NN2 = (mean((Y-Y_fit).^2))^0.5;                            % calculate the RMSE from the optimal NN, for the whole sample (all the points in the simulation, except the discarded ones at the beginning of each thread)
Y_R_squared = 1- ( sum((Y-Y_fit).^2) ) / ( sum((Y-mean(Y)).^2) ); % This computes the R^2 over the whole sample (all the points in the simulation, except the discarded ones at the beginning of each thread)

% ----------------------------------------------
% COMPUTATION OF PLM
% We use NN1_flat to calculate the corresponding PLM - first on fine grid

X1mm=    BB_grid_2D ;
X2mm=    NN_grid_2D ;

X1m = reshape(X1mm,[nval_BB*nval_NN,1]);
X2m = reshape(X2mm,[nval_BB*nval_NN,1]);

X1mr=(X1m-X1r3)/X1r4;
X2mr=(X2m-X2r3)/X2r4;

Xmr = [X1mr X2mr];

Yr_mat2 = f2_NN_eval(Xmr,network_width,NN1_flat);
Y_mat2 = (Yr_mat2*Yr4)+Yr3;
PLM_finegrid_2 = reshape(Y_mat2,size(X1mm));


% We use NN1_flat to calculate the corresponding PLM  - now on HJB grid

X1mm=    squeeze(B(1,1,:,:)) ;
X2mm=    squeeze(N(1,1,:,:)) ;

X1m = reshape(X1mm,[nval_B*nval_N,1]);
X2m = reshape(X2mm,[nval_B*nval_N,1]);

X1mr=(X1m-X1r3)/X1r4;
X2mr=(X2m-X2r3)/X2r4;

Xmr = [X1mr X2mr];

Yr_mat2 = f2_NN_eval(Xmr,network_width,NN1_flat);
Y_mat2 = (Yr_mat2*Yr4)+Yr3;
Y_mat2s= reshape(Y_mat2,size(X1mm));
for ia=1:nval_a
    for iz=1:nval_z
        PLM_2(ia,iz,:,:)=Y_mat2s(:,:);
    end
end

% ALTERNATIVE MACHINE LEARNING SCHEME (ROBUSTNESS): use interpolant to calculate an alternative PLM

PLM_IP = F(X1m,X2m);
if size(Y_fit_IP(:))<10  % On rare ocasions F fails to give a proper output, and I don't want the program to stop because of that
    Y_fit_IP=Y_fit;
    PLM_IP  =PLM_2;
else
    PLM_IP = reshape(PLM_IP,size(X1mm));
end
Y_error_IP = Y-Y_fit_IP;
Y_RMSE_IP = (mean(Y_error_IP.^2))^0.5;


% if we want to use IP instead of NN for the visited area...

% for iB=1:nval_B
%     for iN=1:nval_N
%         for it_sim=1:size(Y,1)
%             if abs(X1(it_sim)-B_grid(iB))<dBB/2
%                 if abs(X2(it_sim)-N_grid(iN))<dNN/2
%                     PLM_2(:,:,iB,iN)=PLM_IP(iB,iN);
%                 end
%             end
%         end
%     end
% end



