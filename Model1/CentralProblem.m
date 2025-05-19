classdef CentralProblem < handle
        % solves the central problem
    properties
        agents
        x_index
        u_index
        nx
        nu
        x0 = [];
        lambda0 = [];
        umin = [];
        umax = [];
        R(1,1) = ones(1);
    end

    methods
        function obj = CentralProblem(agents)
            obj.agents = agents; % {agent_1_central,agent_2_central,agent_3_central};

            obj.x_index = length(obj.agents); % 3 array
            obj.u_index = length(obj.agents); % 3 array

            x_index_count = 0;
            u_index_count = 0;

            for i = 1:length(obj.agents)
                cur_agent = agents{i};
                x_index_count = x_index_count + cur_agent.agentModel.nx_i;  % 1; 1+1; 1+1+1
                u_index_count = u_index_count + cur_agent.agentModel.nu_i;

                obj.x_index(i) = x_index_count; % x_index(1) = Model_1.nx = 1; x_index(2) = 2; x_index(3) = 3
                obj.u_index(i) = u_index_count; % u_index(1) = Model_1.nx = 1; u_index(2) = 2; u_index(3) = 3
   
                % Stack initial conditions
                obj.x0 = [obj.x0;cur_agent.agentModel.x0]; % obj.x0 = [agent_1.x0;agent_2.x0;agent_3.x0] = [0;0;0] 3 x 1
                obj.lambda0 =[obj.lambda0, ones(cur_agent.opt_param.N_d,cur_agent.agentModel.nx_i)]; % 21 x 3
                % stack input constraints
                obj.umin = [obj.umin; cur_agent.agentModel.umin]; % obj.umin = [agent_1.umin;agent_2.umin;agent_3.umin]
                obj.umax = [obj.umax; cur_agent.agentModel.umax]; % 3 x 1

                obj.nx = x_index_count; % 3
                obj.nu = u_index_count; % 3
            end
        end

        function setAgentsStates(obj,x_q,u_q,lambda_q)
            for i=1:length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx_i;
                nui = cur_agent.agentModel.nu_i;
                begin_idx_x = obj.x_index(i); % 1; 2; 3
                end_idx_x = obj.x_index(i) + nxi - 1; % 1; 2; 3

                begin_idx_u = obj.u_index(i); 
                end_idx_u = obj.u_index(i) + nui - 1;

                states.x_i = x_q(begin_idx_x:end_idx_x,:); % x_q = [Model_1.x_q;Model_2.x_q;Model_3.x_q]
                states.u_i = u_q(begin_idx_u:end_idx_u,:); % Extracts the state varibales for the current agent
                states.lambda_i = lambda_q(begin_idx_x:end_idx_x,:);

                cur_agent.setAgentState(states);
            end
        end
%%
        function ui = controlInput(obj,x,lambda,begin_idu,end_idu)
                u0 = (-1/obj.R) * lambda' * x;
                if u0 < obj.umin(begin_idu:end_idu,:)
                    ui = obj.umin(begin_idu:end_idu,:);
                elseif u0 > obj.umax(begin_idu:end_idu,:)
                    ui = obj.umax(begin_idu:end_idu,:);
                else 
                    ui = u0;
                end
        end
%%
         function x_dot = dH_dlambda(obj,t,x,lambda_r,t_span)
            x_dot = zeros(obj.nx,1); % row = 3, col = 1;
            for i = 1 : length(obj.agents) % 1: 3
                cur_agent = obj.agents{i}; 

                begin_idx_x = obj.x_index(i);
                end_idx_x = obj.x_index(i) + cur_agent.agentModel.nx_i - 1;
                begin_idu = obj.u_index(i);
                end_idu = obj.u_index(i) + cur_agent.agentModel.nu_i - 1;   

                lambda_t = interp1(t_span,lambda_r(:,begin_idx_x:end_idx_x),t)'; % 1 x 3
                u_i = obj.controlInput(x(begin_idx_x:end_idx_x,:),lambda_t,begin_idu,end_idu); % x:3 x 1
                x_dot(begin_idx_x:end_idx_x,:) = cur_agent.agentModel.dfdx(x(begin_idx_x:end_idx_x,:),u_i);

                %x_dot(begin_idx_x:end_idx_x) = cur_agent.agentModel.dfdx(x(begin_idx_x:end_idx_x),u_i);
                for j = 1 : length(cur_agent.neighbors)  
                    cur_neighbor = cur_agent.neighbors{j};
                    neighbor_id = cur_neighbor.neighbor_ID;
                    begin_neigh_idx_x = obj.x_index(neighbor_id);
                    end_neigh_idx_x = obj.x_index(neighbor_id) + cur_neighbor.couplingModel.nxj - 1;

                    x_j = x(begin_neigh_idx_x:end_neigh_idx_x,:); % interp1(t_span, x(begin_neigh_idx_x:end_neigh_idx_x,:),t);
                    % select the elements from index begin_idx_x to end_idx_x
                    x_dot(begin_idx_x:end_idx_x,:) = x_dot(begin_idx_x:end_idx_x,:) + cur_neighbor.couplingModel.fij(x(begin_idx_x:end_idx_x,:),x_j);
                end     
            end
         end

          function lambda_dot = dH_dx(obj,t,x,lambda,t_span)
              lambda_dot = zeros(obj.nx,1);

              for i = 1 : length(obj.agents)
                  cur_agent = obj.agents{i};
                    
                  begin_idx_x = obj.x_index(i);
                  end_idx_x = obj.x_index(i) + cur_agent.agentModel.nx_i - 1;

                  begin_idu = obj.u_index(i);
                  end_idu = obj.u_index(i) + cur_agent.agentModel.nu_i - 1;

                  x_i = interp1(t_span,x(:,begin_idx_x:end_idx_x),t);
                  u_i = obj.controlInput(x_i,lambda(begin_idx_x:end_idx_x,:),begin_idu,end_idu);

                  lambda_dot(begin_idx_x:end_idx_x,:) = -cur_agent.agentModel.dldx(x_i) - cur_agent.agentModel.dfdx(x_i,u_i)'*lambda(begin_idx_x:end_idx_x,:);

                  for j = 1 : length(cur_agent.neighbors)
                      cur_neighbor = cur_agent.neighbors{j};
                      neighbor_id = cur_neighbor.neighbor_ID;
                      begin_neigh_idx_x = obj.x_index(neighbor_id);
                      end_neigh_idx_x = obj.x_index(neighbor_id) + cur_neighbor.couplingModel.nxj - 1;
        
                      x_j = interp1(t_span,x(:,begin_neigh_idx_x:end_neigh_idx_x),t);
                      %coupling f_ij l_ij
                      lambda_dot(begin_idx_x:end_idx_x) = lambda_dot(begin_idx_x:end_idx_x,:) - cur_neighbor.couplingModel.dlij_dxi(x_i,x_j)-cur_neighbor.couplingModel.dfij_dxi(x_i,x_j)'*lambda(begin_idx_x:end_idx_x,:);
                      % coupling f_ji l_ji
                      lambda_dot(begin_neigh_idx_x:end_neigh_idx_x) = lambda_dot(begin_neigh_idx_x:end_neigh_idx_x,:) - cur_neighbor.couplingModel.dlij_dxj(x_i,x_j)-cur_neighbor.couplingModel.dfij_dxj(x_i,x_j)'*lambda(begin_idx_x:end_idx_x,:);
                  end
              end    
          end
%%
        function u_q = calcu_uq(obj,x_q,lambda_q)
            [rows, cols] = size(obj.lambda0);
            u_q = zeros(rows, cols);

            for i = 1 : length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx_i;
                nui = cur_agent.agentModel.nu_i;
                begin_idx_x = obj.x_index(i);
                end_idx_x = begin_idx_x + nxi - 1;

                begin_idu_u = obj.u_index(i);
                end_idu_u = begin_idu_u + nui - 1;

                for j = begin_idx_x :end_idx_x
                    for m = 1:size(x_q,1)
                        u_q(m,j) = obj.controlInput(x_q(m,j),lambda_q(m,j),begin_idu_u,end_idu_u);
                    end      
                end    
            end
        end
%%
        function out = dVdx(obj,x_T)
            out= zeros(obj.nx,1);
            for i=1:length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx_i;
                begin_idx_x = obj.x_index(i);
                end_idx_x = obj.x_index(i) + nxi - 1;

                % agent dynamics
                out(begin_idx_x:end_idx_x,:) = cur_agent.agentModel.dVdx(x_T(:,begin_idx_x:end_idx_x));
            end
        end
    end
end