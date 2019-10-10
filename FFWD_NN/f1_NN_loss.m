% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This function evaluates the loss (MSE) for a given set of data

function [y_loss] = f1_NN_loss(y_data,x_data,nwidth,NN,lambda)

nobs   = size(x_data,1);              % number of observations are passed on to here
xwidth = size(x_data,2);              % number of inputs that the NN has, not counting the constant

x_data1 = [ones(nobs,1) x_data];      % add the constant (bias)
xwidth = xwidth+1;

% Construct the NN
w1 = reshape(NN(1:xwidth*nwidth),[xwidth,nwidth]); % size is xwidth x nwidth
w2 = NN(end-nwidth:end-1);                         % size is nwidth x 1
b2 = NN(end);                                      % size is 1x1

% Evaluate the NN for the data and calculate the error
layer_1    = x_data1*w1;             % size is nobs x nwidth
layer_1_NL = log( exp(layer_1)+1 );  % size is nobs x nwidth
y_fitted   = layer_1_NL*w2 + b2;     % size is nobs x 1

% MSE
y_loss     = mean((y_data-y_fitted).^2); % size is 1x1 (this is the MSE)

% Now add the part coming fromt he regularization term

my_regularization = NN;
my_regularization(1:xwidth:nwidth*xwidth)=0;
my_regularization(end)=0;

y_loss     = y_loss + lambda/nobs * (my_regularization'*my_regularization);

end
