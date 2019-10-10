% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This script computes the coefficients of the neural network using line search
% Line search: this algorithm finds a local minimum in the direction of the gradient. This amounts to finding the optimal jumpsize


% STEP 01: COMPUTE THE GRADIENT (this gives us the right direction)
mygrad = f5_NN_gradient(Yrows,Xrows,network_width,NN_flat,lambda); % result is a column vector, with the gradient of the batch

% STEP 1: FIND THREE INITIAL POINTS
% First we need to find three  initial values of jumpsize such that the one in the middle has a lower loss than the ones in the extremes (along the gradient direction)

% Guess 0 : no change
jumpsize0=0;    
NN_flat_0=NN_flat;                                                 % Coefficents of the NN from previous iteration
fval_0=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_0,lambda);     % Evaluate the loss

% Guess 1: jump = learning_speed (from the previous iteration)
jumpsize1=learning_speed_;                                         % First jump is with a predefined step size
NN_flat_1=NN_flat-jumpsize1*mygrad;                                % Gradient descent
fval_1=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_1,lambda);
while fval_1 > fval_0                                              % The first jump should be small enough to ensure a decrease in the loss function 
    jumpsize1=jumpsize1/2;                                         % If the loss function is increasing, reduce the step size as the jump was too big
    NN_flat_1=NN_flat-jumpsize1*mygrad;
    fval_1=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_1,lambda);
end

% Guess 2: Loss in this point should be increasing
jumpsize2=jumpsize1*1.4142;                                        % Now find a point in which loss function is increasing (with respect to guess 1)
NN_flat_2=NN_flat-jumpsize2*mygrad; 
fval_2=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_2,lambda);
while fval_2 < fval_1                                              % If loss function is still decresing, we define this new point as Guess 1 (jumpsize1) and previous Guess 1 as new Guess 0 (jumpsize 0)
    jumpsize0=jumpsize1;
    fval_0 = fval_1;
    NN_flat_0=NN_flat_1;

    jumpsize1=jumpsize2;
    fval_1 = fval_2;
    NN_flat_1=NN_flat_2;

    jumpsize2=jumpsize1*1.1892;                                    % Then we continue looking for a valid Guess 2
    NN_flat_2=NN_flat-jumpsize2*mygrad;
    fval_2=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_2,lambda);
end

% STEP 2: USE A DIVIDE AND CONQUER ALGORITHM TO FIND THE OPTIMUM
% Now we have three points (0,1,2) with the interior one (1) below the other two
% There's at least one local minimum inside there, let's find it by looking
% for points in the middle and discarding the worst ones

for it2 = 1:5
    if fval_2 > fval_0 % Try to get rid of the worst point. Imagine that 2 is worse than 0
        jumpsize8=(jumpsize1+jumpsize2)/2; % mid-point between 1 and 2
        NN_flat_8=NN_flat-jumpsize8*mygrad;
        fval_8=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_8,lambda);
        if fval_8<fval_1       % Either the new point is below 1 -> Replace point 0...
            jumpsize0=jumpsize1;
            fval_0 = fval_1;
            NN_flat_0=NN_flat_1;

            jumpsize1=jumpsize8;
            fval_1 = fval_8;
            NN_flat_1=NN_flat_8;
        else                   % ...Or the new point is above 1 -> Replace point 2
            jumpsize2=jumpsize8;
            fval_2 = fval_8;
            NN_flat_2=NN_flat_8;
        end
    else               % Otherwise 0 is wrose than 2, then we proceed as above but between 0 and 1
        jumpsize8=(jumpsize0+jumpsize1)/2;
        NN_flat_8=NN_flat-jumpsize8*mygrad;
        fval_8=f1_NN_loss(Yrows,Xrows,network_width,NN_flat_8,lambda);
        if fval_8<fval_1
            jumpsize2=jumpsize1;
            fval_2 = fval_1;
            NN_flat_2=NN_flat_1;

            jumpsize1=jumpsize8;
            fval_1 = fval_8;
            NN_flat_1=NN_flat_8;
        else
            jumpsize0=jumpsize8;
            fval_0 = fval_8;
            NN_flat_0=NN_flat_8;
        end
    end
end

% After finishing this, we can say point 1 is a local minimum in the direction of the gradient
