classdef CentralProblem < handle
        % solves the central problem
    properties
        agents
        x_index
        u_index
        nx
        nu
        x0 = [];
        umin = [];
        umax = [];
    end

    methods
        function obj = CentralProblem(agents)
            obj.agents = agents;

            obj.x_index = length(obj.agents);
            obj.u_index = length(obj.agents);

            x_index_count = 1;
            u_index_count = 1;

            for i=1:length(obj.agents)
                cur_agent = agents{i};
                obj.x_index(i) = x_index_count;
                obj.u_index(i) = u_index_count;

                x_index_count = x_index_count + cur_agent.agentModel.nx;
                u_index_count = u_index_count + cur_agent.agentModel.nu;

                % Stack initial conditions
                obj.x0 = [obj.x0;cur_agent.agentModel.x0];

                % stack input constraints
                obj.umin = [obj.umin; cur_agent.agentModel.umin];
                obj.umax = [obj.umax; cur_agent.agentModel.umax];
            end

            obj.nx = x_index_count - 1;
            obj.nu = u_index_count - 1;

        end

        function setAgentsStates(obj,x_q,u_q,lambda_q)
            for i=1:length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx;
                nui = cur_agent.agentModel.nu;
                begin_idx_x = obj.x_index(i);
                end_idx_x = obj.x_index(i) + nxi - 1;

                begin_idx_u = obj.u_index(i);
                end_idx_u = obj.u_index(i) + nui - 1;

                state.x_i = x_q(begin_idx_x:end_idx_x,:);
                state.u_i = u_q(begin_idx_u:end_idx_u,:);
                state.lambda_i = lambda_q(begin_idx_x:end_idx_x,:);

                cur_agent.setAgentState(state);
            end


        end

        function x_dot = dH_dlambda(obj,t,x,u,t_span)
            u = interp1(t_span,u',t)';
            x_dot = zeros(obj.nx,1);

            % iterate over agents
            for i=1:length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx;
                nui = cur_agent.agentModel.nu;
                begin_idx_x = obj.x_index(i);
                end_idx_x = obj.x_index(i) + nxi - 1;

                begin_idx_u = obj.u_index(i);
                end_idx_u = obj.u_index(i) + nui - 1;

                % agent dynamics
                x_dot(begin_idx_x:end_idx_x) = x_dot(begin_idx_x:end_idx_x) +...
                    cur_agent.agentModel.ffct(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u));

                % iterate over neighbors
                for j=1:length(cur_agent.neighbors)
                    cur_neighbor = cur_agent.neighbors{j};
                    nxj = cur_neighbor.couplingModel.nxj;
                    nuj = cur_neighbor.couplingModel.nuj;
                    neighbor_id = cur_neighbor.neighbor_ID;

                    begin_neighbor_idx_x = obj.x_index(neighbor_id);
                    end_neighbor_idx_x = obj.x_index(neighbor_id) + nxj - 1;

                    begin_neighbor_idx_u = obj.u_index(neighbor_id);
                    end_neighbor_idx_u = obj.u_index(neighbor_id) + nuj - 1;

                    %coupling dynamics
                    x_dot(begin_idx_x:end_idx_x) = x_dot(begin_idx_x:end_idx_x) + ...
                        cur_neighbor.couplingModel.fij(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                        x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u));
                end

            end
        end


        function lambda_dot = dH_dx(obj,t,x,u,lambda,t_span)
            x = interp1(t_span,x',t)';
            u = interp1(t_span,u',t)';
            lambda_dot = zeros(obj.nx,1);

            % iterate over agents
            for i=1:length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx;
                nui = cur_agent.agentModel.nu;

                begin_idx_x = obj.x_index(i);
                end_idx_x = obj.x_index(i) + nxi - 1;

                begin_idx_u = obj.u_index(i);
                end_idx_u = obj.u_index(i) + nui - 1;

                % local adjoint dynamics
                lambda_dot(begin_idx_x:end_idx_x) =  lambda_dot(begin_idx_x:end_idx_x) - (...
                    cur_agent.agentModel.dldx(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u)) + ...
                    cur_agent.agentModel.dfdx(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u))'*lambda(begin_idx_x:end_idx_x));

                % iterate over neighbors
                for j=1:length(cur_agent.neighbors)
                    cur_neighbor = cur_agent.neighbors{j};
                    nxj = cur_neighbor.couplingModel.nxj;
                    nuj = cur_neighbor.couplingModel.nuj;
                    neighbor_id = cur_neighbor.neighbor_ID;

                    begin_neighbor_idx_x = obj.x_index(neighbor_id);
                    end_neighbor_idx_x = obj.x_index(neighbor_id) + nxj - 1;

                    begin_neighbor_idx_u = obj.u_index(neighbor_id);
                    end_neighbor_idx_u = obj.u_index(neighbor_id) + nuj - 1;

                    %coupling dynamics
                    lambda_dot(begin_idx_x:end_idx_x) =  lambda_dot(begin_idx_x:end_idx_x) - ( ...
                        cur_neighbor.couplingModel.dlij_dxi(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                        x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u)) + ...
                        cur_neighbor.couplingModel.dfij_dxi(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                        x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u))' * lambda(begin_idx_x:end_idx_x));

                    %coupling dynamics
                    lambda_dot(begin_neighbor_idx_x:end_neighbor_idx_x) = lambda_dot(begin_neighbor_idx_x:end_neighbor_idx_x) - ( ...
                        cur_neighbor.couplingModel.dlij_dxj(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                        x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u)) + ...
                        cur_neighbor.couplingModel.dfij_dxj(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                        x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u))' * lambda(begin_idx_x:end_idx_x));
                 end

            end
        end

        function grad_u = dH_du(obj,t,x,u,lambda)
            grad_u = zeros(obj.nu,length(t));

            % Iterate over horizon
            for k=1:length(t)
                x_n = x(:,k);
                u_n = u(:,k);
                lambda_n = lambda(:,k);
                % iterate over agents
                for i=1:length(obj.agents)
                    cur_agent = obj.agents{i};
                    nxi = cur_agent.agentModel.nx;
                    nui = cur_agent.agentModel.nu;

                    begin_idx_x = obj.x_index(i);
                    end_idx_x = obj.x_index(i) + nxi - 1;

                    begin_idx_u = obj.u_index(i);
                    end_idx_u = obj.u_index(i) + nui - 1;

                    % local gradient
                    grad_u(begin_idx_u:end_idx_u,k) = grad_u(begin_idx_u:end_idx_u,k) + (...
                        cur_agent.agentModel.dldu(t,x_n(begin_idx_x:end_idx_x),u_n(begin_idx_u:end_idx_u)) + ...
                        cur_agent.agentModel.dfdu(t,x_n(begin_idx_x:end_idx_x),u_n(begin_idx_u:end_idx_u))'*lambda_n(begin_idx_x:end_idx_x));

                    % iterate over neighbors
                    for j=1:length(cur_agent.neighbors)
                        cur_neighbor = cur_agent.neighbors{j};
                        nxj = cur_neighbor.couplingModel.nxj;
                        nuj = cur_neighbor.couplingModel.nuj;
                        neighbor_id = cur_neighbor.neighbor_ID;

                        begin_neighbor_idx_x = obj.x_index(neighbor_id);
                        end_neighbor_idx_x = obj.x_index(neighbor_id) + nxj - 1;

                        begin_neighbor_idx_u = obj.u_index(neighbor_id);
                        end_neighbor_idx_u = obj.u_index(neighbor_id) + nuj - 1;


                        %coupling dynamics
                        grad_u(begin_idx_u:end_idx_u,k) =  grad_u(begin_idx_u:end_idx_u,k) + ( ...
                            cur_neighbor.couplingModel.dlij_dui(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                            x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u)) + ...
                            cur_neighbor.couplingModel.dfij_dui(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                            x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u))' * lambda_n(begin_idx_x:end_idx_x));

                        %coupling dynamics
                        grad_u(begin_neighbor_idx_u:end_neighbor_idx_u,k) = grad_u(begin_neighbor_idx_u:end_neighbor_idx_u,k) + ( ...
                            cur_neighbor.couplingModel.dlij_duj(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                            x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u))+ ...
                            cur_neighbor.couplingModel.dfij_duj(t,x(begin_idx_x:end_idx_x),u(begin_idx_u:end_idx_u),...
                            x(begin_neighbor_idx_x:end_neighbor_idx_x),u(begin_neighbor_idx_u:end_neighbor_idx_u))' * lambda_n(begin_idx_x:end_idx_x));
                    end
                end
            end
        end

        function out = dVdx(obj,T,x_T)

            out= zeros(obj.nx,1);
            for i=1:length(obj.agents)
                cur_agent = obj.agents{i};
                nxi = cur_agent.agentModel.nx;

                begin_idx_x = obj.x_index(i);
                end_idx_x = obj.x_index(i) + nxi - 1;

                % agent dynamics
                out(begin_idx_x:end_idx_x) = cur_agent.agentModel.dVdx(T,x_T(begin_idx_x:end_idx_x));

            end
        end

    end
end