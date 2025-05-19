classdef couplingModel_2_1 < handle
    % Dynamics of a double integrator

    properties
        nxi{mustBeInteger} = 1;                                   % Number of states      
        nxj{mustBeInteger} = 1;                                   % Number of states of neighbor
        modelparams(:,1);
        P(1,1) = ones(1);
    end

    methods
        function obj = couplingModel_2_1(params)
            obj.modelparams = params.modelparams;
        end

        % System dynamics
        function x_dot = fij(obj,x_i,x_j)
            %eps = obj.modelparams(1);
            x_dot = cos(x_j - x_i);
        end

        % Jacobian w.r.t.the agent states of system dynamics
        function jac_x_i = dfij_dxi(obj,x_i,x_j)
            %eps = obj.modelparams(1);
            jac_x_i = sin(x_j - x_i);
        end

        function jac_x_j = dfij_dxj(obj,x_i,x_j)
            %eps = obj.modelparams(1);
            jac_x_j = -sin(x_j - x_i);
        end


        % cost function
        function cost = lijfct(obj,x_i,x_j)
            cost = 0;
        end

        % Gradient w.r.t. the states of cost function
        function grad_x_i = dlij_dxi(obj,x_i,x_j)
            grad_x_i = (0);
        end

        function grad_x_j = dlij_dxj(obj,x_i,x_j)
            grad_x_j = (0);
        end
    end
end