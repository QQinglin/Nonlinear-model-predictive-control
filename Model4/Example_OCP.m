clear;

% Optimization Parameters 
opt_param.N_d = 21;                                    % Discretization Points
opt_param.T = 1;                                      % Prediction horizon
opt_param.t = linspace(0,opt_param.T,opt_param.N_d);  % Uniform time grid                          

% Sensi Parameters 
opt_param.sens_iter_k = 30;  
opt_param.fixed_iter_j = 30;

% Agent Model parameters
agentModel_param.P = diag((1));
agentModel_param.Q = diag((1));
agentModel_param.R = 1;                           
agentModel_param.umin = -5;
agentModel_param.umax = 5;

% Coupling Model parameters
couplingModel_param.modelparams = 1.0;

% create models and coupling models object
agent_model_1 = agentModel_1(agentModel_param);
agent_model_2 = agentModel_2(agentModel_param);
agent_model_3 = agentModel_3(agentModel_param);
coupling_model_1_2 = couplingModel_1_2(couplingModel_param);
coupling_model_2_1 = couplingModel_2_1(couplingModel_param);
coupling_model_2_3 = couplingModel_2_3(couplingModel_param);
coupling_model_3_2 = couplingModel_3_2(couplingModel_param);

% Define Coupling Graph, adjecency matrix with couplings
opt_param.A = {[],coupling_model_1_2,[];coupling_model_2_1,[],coupling_model_2_3;[],coupling_model_3_2,[]};  

% Create Agents object
agent_id = 1; 
agentModel_1 = Agent(opt_param,agent_model_1,agent_id); 
agent_id = 2; 
agentModel_2 = Agent(opt_param,agent_model_2,agent_id);  
agent_id = 3;
agentModel_3 = Agent(opt_param,agent_model_3,agent_id);  
agents = {agentModel_1,agentModel_2,agentModel_3};

%% Solve Central OCP
%Central model parameters
model_param.P = diag([1,1,1]);
model_param.Q = diag([1,1,1]);
model_param.R = diag([1,1,1]);
model_param.x0 = [1;1;1]; 
model_param.umin = -5;
model_param.umax = 5;
%set Model
model = central_model(model_param);
% Optimization Parameters 
% opt_param.N_d = 21;                             
opt_param.lambda0 = ones(opt_param.N_d,model.nx);  
[u_k,lambda_k,x_k] = central_solve(opt_param, model);

%%%%%%%%%%%%%%%%%%%
% run Algorithm 
%%%%%%%%%%%%%%%%%%%
 run_sensi(opt_param,agents,u_k,lambda_k,x_k);


