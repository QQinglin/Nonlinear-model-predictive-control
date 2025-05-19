classdef Agent < handle
    % Takes care of all the local steps
    properties
        agentModel;
        x_q;
        u_q;
        lambda_q;
        agent_ID;
        neighbors;
        opt_param;

        last_x = [];
        first_x = [];
        last_lambda = [];
        collect_result = [];

        % Anderson acceleration
        conver_k = [];
        sen_k = [];
        last_re = [];



       
        x = [];
        x_bar = [];
    end
    methods
        % Constructor function for the class agent initializes obj.
        function obj = Agent(param,model,global_id)
            % global Id for identification
            obj.agent_ID = global_id;                 

            % agent Model
            obj.agentModel = model;                   

            % Sensi Initialization
            obj.x_q = zeros(param.N_d,obj.agentModel.nx_i);      
            obj.lambda_q = zeros(param.N_d,obj.agentModel.nx_i);  
            obj.u_q = zeros(param.N_d,obj.agentModel.nu_i);    

            obj.x = ones(param.N_d,obj.agentModel.nx_i);
            obj.last_re = ones(param.N_d,obj.agentModel.nx_i);

            % opt params
            obj.opt_param = param;                    

            % Initialize neighbors
            % couplings via adjacency matrix
            neighbors_idx = find(~cellfun(@isempty,param.A(global_id,:)));   % array [1 3] index of each neighbor
            num_neighbors = length(neighbors_idx);  % 2
            
            % creates a cell with empty element ,num_neighbors specifies the number of rows, and 1 specifies the number of columns. 
            obj.neighbors = cell(num_neighbors,1);        % {[];[]}  
            % Create neighbors
            for j=1:num_neighbors   
               % loop through each index of neighbors，and everytime creating an object of class <Neighbor>.
               % param.A{global_id,neighbors_idx(j)} specifies the index of coupling models in matrix A
               % param initinize object x_j, lambda_j [0 0 0...]
               % neighbors_idx(j) specifies the neighbor_ID
               % the object is assigned to the j-th element of neighbors {[],[]}
               % the cell arrray looks like this {[obj.neigh_1],[obj.neigh_3]},first_neighbor = obj.neighbors{1};
               % use such as first_neighbor_x_j = obj.neighbors{1}.x_j;commands to access and manipulate the properties of obj
                obj.neighbors{j} = Neighbor(param,param.A{global_id,neighbors_idx(j)},neighbors_idx(j)); 
            end
        end

        function compute_sensitivities(obj)
            x_i = obj.x_q;                       
            
            % Compute sensitivities for all neighbors
            for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};   
                sensi_coupling = obj.opt_param.A{cur_neighbor.neighbor_ID,obj.agent_ID};
                x_j = cur_neighbor.x_j;           
                lambda_j = cur_neighbor.lambda_j;  
                
                cur_neighbor.gji_x = sensi_coupling.dlij_dxj(x_j,x_i) ...
                        +  sensi_coupling.dfij_dxj(x_j,x_i) .* lambda_j;
            end
        end
%%
        function states = getAgentState(obj)   % get state from OCP solution
                 states = struct();
                 states.x_i = obj.x_q;  
                 states.lambda_i = obj.lambda_q; 
        end

        function setNeighborsStates(obj,states,from_Id) 
            % send trajectories to all neighbors
             for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};   
                if(cur_neighbor.neighbor_ID == from_Id)
                   cur_neighbor.x_j = states.x_i; 
                   cur_neighbor.lambda_j = states.lambda_i;
                    break
                end 
             end                  
        end
%%
        function x_dot = dHi_dlambda_i(obj,t,x_i,lambda_r,t_span)
            lambda_t = interp1(t_span,lambda_r,t)';
            u_i = obj.agentModel.controlInput(x_i,lambda_t);

            x_dot = obj.agentModel.ffct(x_i,u_i);

            % influence of neighbor
            for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};
                x_j = interp1(t_span,cur_neighbor.x_j',t)';
                x_dot = x_dot + cur_neighbor.couplingModel.fij(x_i,x_j);   % fij(obj,~,x_i,~,x_j,~) = x_j(1) - x_i(2);
            end
        end

        function lambda_dot = dHi_dxi(obj,t,x_r,lambda_i,t_span)
            x_i = interp1(t_span,x_r',t)';
            u_i = obj.agentModel.controlInput(x_i,lambda_i);

            lambda_dot = -(obj.agentModel.dldx(x_i) + obj.agentModel.dfdx(x_i,u_i)' * lambda_i);

            % influence of neighbor
            for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};
                x_j = interp1(t_span,cur_neighbor.x_j,t);
                gji_x = interp1(t_span,cur_neighbor.gji_x,t); 

                % part due to coupled dynamics
                lambda_dot = lambda_dot - (cur_neighbor.couplingModel.dlij_dxi(x_i,x_j) + ...
                    (cur_neighbor.couplingModel.dfij_dxi(x_i,x_j)' * lambda_i));

                % part due to sensitivities
                lambda_dot = lambda_dot - gji_x;
            end
        end

%%
        function solve_localOCP(obj)
            % solves the OCP via fixed point iteration method
            x0 = obj.agentModel.x0;
            lambda_r = obj.lambda_q;
            t_span = obj.opt_param.t;
            T_end = t_span(end);
            %% Begin of local OCP solution via fixed point iteration method 
            for r = 1:obj.opt_param.fixed_iter_j
                % Forward integration of the state dynamics
                x_r = ode4(@(t,x)obj.dHi_dlambda_i(t,x,lambda_r,t_span),t_span,x0);

                if r == 1
                    obj.first_x = [obj.first_x,x_r];
                end
                % backward integration of the adjoint state dynamics
                rev_t_span = flip(t_span);
                x_T = x_r(end, :)';
                lambda_T = obj.agentModel.dVdx(x_T);
                lambda_r_rev = ode4(@(t,lambda)obj.dHi_dxi(t,x_r,lambda,t_span),rev_t_span,lambda_T);
                lambda_r = flip(lambda_r_rev);
            end

            obj.last_lambda = [obj.last_lambda,lambda_r];

            obj.last_x = [obj.last_x,x_r];
    
            % update the output
            u_r = zeros(size(obj.u_q,1),1);
            for i = 1 : size(obj.u_q,1)
                x_current = x_r(i,:)';
                lambda_current = lambda_r(i,:)';
                u_r(i) = obj.agentModel.controlInput(x_current,lambda_current);
            end  
             % update the state of object
            obj.x_q = x_r;  
            obj.lambda_q = lambda_r;
            obj.u_q = u_r;
        end
    end
end

