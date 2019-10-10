% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This function evaluates the NN for a given set of data

function [y_fitted] = f2_NN_eval(x_data,nwidth,NN)

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

end
