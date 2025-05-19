clear;

% Optimization Parameters 
opt_param.N_d = 41;                                   % Discretization Points
opt_param.T = 1;                                    % Prediction horizon
opt_param.t = linspace(0,opt_param.T,opt_param.N_d);  % Uniform time grid
opt_param.grad_iter = 10;                             % Maximum number of gradient iterations

% Sensi Parameters 
opt_param.sensi_iter = 30;

% Agent Model parameters
agentModel_param.P = diag([10]);
agentModel_param.Q = diag([10]);
agentModel_param.R = 1;
agentModel_param.x0 = [0.5];                            
agentModel_param.umin = -1;
agentModel_param.umax = 1;

% Coupling Model parameters
couplingModel_param.modelparams = 1.0;

% set agent and coupling models

agent_model = DoubleIntegrator_agentModel(agentModel_param);

% x_dot = [x_1_dot;x_2_dot;x_3_dot];
coupling_model = DoubleIntegrator_couplingModel(couplingModel_param);

% Sensi Initialization 
init.x_0 = ones(agent_model.nx,opt_param.N_d); 
init.lambda_0 = ones(agent_model.nx,opt_param.N_d);

% Create Agents
% 每一次创建Agent ，传入的模型都是相同的，并且是 整个系统的动态方程而不是子系统的
agent_id = 1;
dynamics_func = @(x_i, u_i) (x_i(1) + u_i(1) * x_i(1));
dynamics_jac = @(x_i, u_i) [1+u(1),0,0];
agent_model_1 = DoubleIntegrator_agentModel(agentModel_param,dynamics_func,dynamics_jac);
coupl_fun = @(x_i, x_j) (x_j(2) - x_i(1));
coupl_jac = @(x_i, u_i) [-1,0,0];
coupl_cost = @(x_i, u_i) (0);
coupl_cost_jac  = @(x_i, u_i) [0,0,0];
coupling_model_12 = DoubleIntegrator_couplingModel(couplingModel_param,coupl_fun,coupl_jac,coupl_cost,coupl_cost_jac);

agent_id = 2;
dynamics_func = @(x_i, u_i) (2 * x_i(2) + u_i(2) * x_i(2));
dynamics_jac = @(x_i, u_i) [0,2+u(2),0];
agent_model_2 = DoubleIntegrator_agentModel(agentModel_param,dynamics_func,dynamics_jac);
coupl_fun = @(x_i, x_j) (x_j(1) - x_i(2) + x_j(3) - x_i(2));
coupl_jac = @(x_i, u_i) [0,-1,0];
coupl_cost = @(x_i, u_i) (0);
coupl_cost_jac  = @(x_i, u_i) [0,0,0];
coupling_model_21 = DoubleIntegrator_couplingModel(couplingModel_param,coupl_fun,coupl_jac,coupl_cost,coupl_cost_jac);
coupl_fun = @(x_i, x_j) (x_j(3) - x_i(2));
coupl_jac = @(x_i, u_i) [0,-1,0];
coupl_cost = @(x_i, u_i) (0);
coupl_cost_jac  = @(x_i, u_i) [0,0,0];
coupling_model_23 = DoubleIntegrator_couplingModel(couplingModel_param,coupl_fun,coupl_jac,coupl_cost,coupl_cost_jac);

agent_id = 3;
dynamics_func = @(x_i, u_i) (3 * x_i(3) + u_i(3) * x_i(3));
dynamics_jac = @(x_i, u_i) [0,0,3+u(3)];
agent_model_3 = DoubleIntegrator_agentModel(agentModel_param,dynamics_func,dynamics_jac);
coupl_fun = @(x_i, x_j) (x_j(2) - x_i(3));
coupl_jac = @(x_i, u_i) [0,0,-1];
coupl_cost = @(x_i, u_i) (0);
coupl_cost_jac  = @(x_i, u_i) [0,0,0];
coupling_model_32 = DoubleIntegrator_couplingModel(couplingModel_param,coupl_fun,coupl_jac,coupl_cost,coupl_cost_jac);

opt_param.A = {[],coupling_model_12,[];coupling_model_21,[],coupling_model_23;[],coupling_model_32,[]};  % Adjecency matrix with couplings

agent_1_sensi = Agent(init,opt_param,agent_model_1,agent_id);
agent_2_sensi = Agent(init,opt_param,agent_model_2,agent_id);
agent_3_sensi = Agent(init,opt_param,agent_model_3,agent_id);

agents = {agent_1_sensi,agent_2_sensi,agent_3_sensi};



%%%%%%%%%%%%%%%%%%%
% run Algorithm (Assumption f_ij = f_ji)!!!!!!!
%%%%%%%%%%%%%%%%%%
run_sensi(opt_param,agents);

