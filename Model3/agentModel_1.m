classdef agentModel_1 < handle

    properties
        nx_i{mustBeInteger} = 1;               % Number of states
        nu_i{mustBeInteger} = 1;               % Number of controls
        P(1,1) = ones(1);
        Q(1,1) = ones(1);
        R(1,1) = ones(1);
        umin(1,1) = -inf;
        umax(1,1) = inf;
        x0(1,1)   = ones(1,1);
        %modelparams(:,1);
    end

    methods
        function obj = agentModel_1(params)

            obj.P = params.P;
            obj.Q = params.Q;
            obj.R = params.R;
            obj.umin = params.umin;
            obj.umax = params.umax;
        end

        % System dynamics
        function x_dot = ffct(obj,x_i,u_i)
            x_dot = -x_i^3 + u_i*x_i;
        end

        % Jacobian w.r.t.the states of system dynamics
        function jac_x = dfdx(obj,x_i,u_i)
             jac_x = [-3*x_i^2 + u_i];
        end

        function ui = controlInput(obj,x,lambda)
            u0 = (-1/obj.R) * lambda' * x;
            if u0 < obj.umin
                ui = obj.umin;
            elseif u0 > obj.umax
                ui = obj.umax;
            else 
                ui = u0;
            end
        end

        % cost function
        function cost = lfct(obj,x_i,u_i)
            cost = 0.5 * x_i' * obj.Q * x_i + 0.5 * u_i'*obj.R * u_i;
        end

        % Gradient w.r.t. the states of cost function
        function grad_x = dldx(obj,x_i)
            grad_x = obj.Q * x_i;
        end

        % terminal cost
        function cost = Vfct(obj,x_T)
            cost = 0.5 * x_T' * obj.P * x_T;
        end

        % Gradient w.r.t. to the states of terminal cost
        function grad_x = dVdx(obj,x_T)
            grad_x = obj.P*x_T;
        end
    end
end