classdef couplingModel_2_1 < handle
    % Dynamics of a double integrator

    properties
        nxi{mustBeInteger} = 3;                                   % Number of states      
        nxj{mustBeInteger} = 3;                                   % Number of states of neighbor
        modelparams(:,1);
        P(1,1) = ones(1);
        x_ijref;
        Q;
    end

    methods
        function obj = couplingModel_2_1(params)
            obj.modelparams = params.modelparams;
            obj.Q = diag(ones(1,3));
            obj.x_ijref = [1,0,0]';   % assume keeping the distance delta_x = -1, agent 2 is behind agent 1
            % obj.x_ijref = [ones(21, 1), zeros(21, 2)];
        end

        % System dynamics
        function x_dot = fij(obj,x_i,x_j)
            %eps = obj.modelparams(1);
            x_dot = zeros(3,1);
        end

        % Jacobian w.r.t.the agent states of system dynamics
        function jac_x_i = dfij_dxi(obj,x_i,x_j)
            %eps = obj.modelparams(1);
            jac_x_i = zeros(3,3);
        end

        function jac_x_j = dfij_dxj(obj,x_i,x_j)
            %eps = obj.modelparams(1);
            jac_x_j = zeros(3,3);   % **********************************
        end

        % cost function
        function cost = lijfct(obj,x_i,x_j) % x_i is varible of agent 2
            cost = 0.5 * (x_j - x_i + obj.x_ijref)' * obj.Q * (x_j - x_i + obj.x_ijref);
        end

        % Gradient w.r.t. the states of cost function
        function grad_x_i = dlij_dxi(obj,x_i,x_j)
            grad_x_i = -obj.Q * (x_j - x_i + obj.x_ijref);
        end

        function grad_x_j = dlij_dxj(obj,x_i,x_j)
            grad_x_j = obj.Q * (x_j' - x_i' + obj.x_ijref); % Q 3x3 unite matrix, x_j - x_i + obj.x_ijref 21x3
        end % grad_x_j 3x21
    end
end