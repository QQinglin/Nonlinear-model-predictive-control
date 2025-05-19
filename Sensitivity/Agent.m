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
        conv_u_q;
    end

    methods
        function obj = Agent(init,param,model,global_id)

            % global Id for identification
            obj.agent_ID = global_id;

            % init_quantities
            obj.x_q = init.x_0;
         
            obj.lambda_q = init.lambda_0;

            % agent Model
            obj.agentModel = model;

            % opt params
            obj.opt_param = param;


            % Initialize neighbors
            % couplings via adjacency matrix
            neighbors_idx = find(~cellfun(@isempty,param.A(global_id,:)));   % the second row of A ; @isempty get [true, false, true];find get [1, 3]

            num_neighbors = length(neighbors_idx); % how many neighbors, here is 2
            obj.neighbors = cell(num_neighbors,1); % create a array, here is {[], []}

            % Create neighbors
            for j=1:num_neighbors
                % j=1 --> obj.neighbors{1} = Neighbor(param, param.A{2, 1}, 1);
                % j=2 --> obj.neighbors{2} = Neighbor(param, param.A{2, 3}, 3);
                obj.neighbors{j} = Neighbor(param,param.A{global_id,neighbors_idx(j)},neighbors_idx(j));
            end
        end
%%
        function states = getAgentState(obj)     % cur_agent.getAgentState() , cur_agent = 2
                 states.x_i = obj.x_q;
                 %states.u_i = obj.u_q;
                 states.lambda_i = obj.lambda_q;
        end

        function setAgentState(obj,states)
                 obj.x_q = states.x_i;
                 obj.u_q = states.u_i;
                 obj.lambda_q = states.lambda_i;
        end

        function setNeighborsStates(obj,states,cur_agent_Id)  % agents{1}.setNeighborsStates(cur_agent.getAgentState(),2);  cur_agent = 2
             for j = 1:length(obj.neighbors)               
                cur_neighbor = obj.neighbors{j};   % obj.neighbors{2} = Neighbor(param, param.A{2, 3}, 1 or 3);    param.A {global_id,neighbors_idx(j)}
                if(cur_neighbor.neighbor_ID == cur_agent_Id)
                   cur_neighbor.x_j = states.x_i;
                   cur_neighbor.lambda_j = states.lambda_i;
                    break
                end 
             end                  
        end

        function x_dot = dHi_dlambdai(obj,t, x, lambda_r,t_span)
            lambda_t = interp1(t_span, lambda_r, t);
            u_i = x * lambda_t;
            x_dot = obj.agentModel.ffct(x,u_i);

            % influence of neighbor
            for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};
                x_j = interp1(t_span,cur_neighbor.x_j',t)';
                x_dot = x_dot + cur_neighbor.couplingModel.fij(x,x_j);
            end
        end

        function lambda_dot = dHi_dxi(obj,t,x_r,lambda,t_span)
            x_i = interp1(t_span,x_r',t)';
            u_i = x_i * lambda;

            % local part
            lambda_dot = obj.agentModel.dldx(t,x_i,u_i) + obj.agentModel.dfdx(t,x_i,u_i)' * lambda_i;

            % influence of neighbor
            for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};
                x_j = interp1(t_span,cur_neighbor.x_j',t)';
                gji_x = interp1(t_span,cur_neighbor.gji_x',t)';

                % part due to coupled dynamics
                lambda_dot = -lambda_dot - (cur_neighbor.couplingModel.dlij_dxi(x_i,x_j) + ...
                cur_neighbor.couplingModel.dfij_dxi(x_i,x_j)' * lambda_i);
            end
        end
%%
        function compute_sensitivities(obj)
            x_i = obj.x_q;

            % Compute sensitivities for all neighbors
            for j = 1:length(obj.neighbors)
                cur_neighbor = obj.neighbors{j};
                x_j = cur_neighbor.x_j;
                lambda_j = cur_neighbor.lambda_j;

                % pointwise in time
                for k = 1:obj.opt_param.N_d
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % (Assumption f_ij = f_ji, l_ij = l_ji)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % gradient w.r.t. x_i of neighbors cost functional
                    cur_neighbor.gji_x(:,k) = cur_neighbor.couplingModel.dlij_dxi(x_i(:,k),x_j(:,k)) ...
                        +  cur_neighbor.couplingModel.dfij_dxi(x_i(:,k),x_j(:,k))' * lambda_j(:,k);
                end

            end
        end
%% Fixed-point-iteration beginns
 function solve_localOCP(obj)
            % warm start with initial input from last iteration
            % initialstate, time distribution, step size
            x0 = ones(agent_model.nx,1);
            lambda_r = obj.lambda_q;
            t_span = obj.opt_param.t;
            T_end = t_span(end);

            for r=1:obj.opt_param.grad_iter
                x_r = ode4(@(t,x)obj.dHi_dlambdai(t,x,lambda_r,t_span),t_span,x0)';
                % backward integration of the adjoint state dynamics
                rev_t_span = flip(t_span);
                x_T = x_r(:,end);
                lambda_T = obj.agentModel.dVdx(T_end,x_T);
                lambda_r_rev = ode4(@(t,lambda)obj.dHi_dxi(t,x_r,lambda,t_span),rev_t_span,lambda_T);
                lambda_r = flip(lambda_r_rev)';
            end

            % output
        end

    end
end

