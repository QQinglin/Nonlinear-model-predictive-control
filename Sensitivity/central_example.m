clear;close all;

% Model parameters
model_param.P = diag([1,1,1]);
model_param.Q = diag([1,1,1]);
model_param.R = diag([1,1,1]);
model_param.x0 = [1;1;1]; 
model_param.umin = -5;
model_param.umax = 5;
%set Model
model = Model(model_param);

% Optimization Parameters 
opt_param.N_d = 21;                                 % Discretization Points
opt_param.T = 1;                                   % Prediction horizon
opt_param.t = linspace(0,opt_param.T,opt_param.N_d); % Uniform time grid
opt_param.lambda0 = ones(opt_param.N_d,model.nx);  

[u_k,lambda_k,x_k] = solve(opt_param, model);
