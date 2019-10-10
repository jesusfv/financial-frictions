% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This function computes the gradient using the backpropagation algorithm

function [my_gradient] = f5_NN_gradient(y_data,x_data,nwidth,NN,lambda)

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

y_error    = y_fitted-y_data;        % size is nobs x 1

% Third, back-prop
y_error_wide = repmat(y_error,[1,nwidth]);  % size is nobs x nwidth
w2_grad      = y_error_wide.*layer_1_NL;    % size is nobs x nwidth
b2_grad      = y_error;                     % size is nobs x 1


layer_1_grad = 1./(1+exp(-layer_1));        % size is nobs x nwidth

w1_grad=zeros(nobs,xwidth,nwidth);
for it0 = 1:nobs
    for it1=1:xwidth
        for it2=1:nwidth
            w1_grad(it0,it1,it2) = layer_1_grad(it0,it2)*w2(it2,1)*y_error(it0,1)*x_data1(it0,it1);
        end
    end
end

w1_grad_mean = mean(w1_grad,1);
w2_grad_mean = mean(w2_grad,1);
b2_grad_mean = mean(b2_grad,1);

% Finally, put gradient into tensor form
my_gradient = [ w1_grad_mean(:) ; w2_grad_mean(:) ; b2_grad_mean(:) ];

% Regularization added for every weight (but not for the bias coefficients, regularization is not usually applied there)

my_regularization = NN;
my_regularization(1:xwidth:nwidth*xwidth)=0;
my_regularization(end)=0;

my_gradient = 2*( my_gradient + (lambda/nobs)*my_regularization );


end
