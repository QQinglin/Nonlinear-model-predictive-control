clear;

% Optimization Parameters 
opt_param.N_d = 21;                                    % Discretization Points
opt_param.T = 2;                                      % Prediction horizon
opt_param.t = linspace(0,opt_param.T,opt_param.N_d);  % Uniform time grid                          

% Sensi Parameters 
opt_param.sens_iter_k = 30;  
opt_param.fixed_iter_j = 30;

% Agent Model parameters
agentModel_param.P = diag(ones(1, 3));
agentModel_param.Q = diag(ones(1, 3));
agentModel_param.R = diag(ones(1, 2));                          
agentModel_param.umin = [-5;-5];
agentModel_param.umax = [5;5];

% Coupling Model parameters
couplingModel_param.modelparams = 1.0;

% create models and coupling models object
agent_model_1 = agentModel_1(agentModel_param);
agent_model_2 = agentModel_2(agentModel_param);

coupling_model_1_2 = couplingModel_1_2(couplingModel_param);
coupling_model_2_1 = couplingModel_2_1(couplingModel_param);

% Define Coupling Graph, adjecency matrix with couplings
opt_param.A = {[],coupling_model_1_2;
               coupling_model_2_1,[]};  

% Create Agents object
agent_id = 1; 
agentModel_1 = Agent(opt_param,agent_model_1,agent_id); 
agent_id = 2; 
agentModel_2 = Agent(opt_param,agent_model_2,agent_id);  
agents = {agentModel_1,agentModel_2};

%%%%%%%%%%%%%%%%%%%
% run Algorithm 
%%%%%%%%%%%%%%%%%%
 % run_sensi(opt_param,agents);

%% Solve Central OCP
%Central model parameters
model_param.P = diag([1,1,1]);
model_param.Q = diag([1,1,1]);
model_param.R = diag([1,1]);
model_param.x0 = [0;0;0]; 
model_param.umin = [-1;-1];
model_param.umax = [3;3];
%set Model
model = central_model(model_param);
% Optimization Parameters 
% opt_param.N_d = 21;  
model_param.t = opt_param.t;
model_param.fixed_iter_j = opt_param.fixed_iter_j;
model_param.lambda0 = ones(opt_param.N_d,model.nx);  
solve_centralOCP(model_param, model); % [u_k,lambda_k,x_k] = solve_centralOCP(opt_param, model);
