classdef agentModel_2 < handle

    properties
        nx_i{mustBeInteger} = 3;               % Number of states
        nu_i{mustBeInteger} = 2;               % Number of controls
        P = diag(ones(1, 3));
        Q = diag(ones(1, 3));
        R = diag(ones(1, 2));
        umin = [];
        umax = [];
        x0 = [1,1,1]'; % initial conditions: p_x = 0, p_y = 1, yaw angle = 1
        %modelparams(:,1);

        x_iref;
    end

    methods
        function obj = agentModel_2(params)

            obj.P = params.P;
            obj.Q = params.Q;
            obj.R = params.R;
            obj.umin = params.umin;
            obj.umax = params.umax;

            obj.x_iref = [4,4,5]'; % set p_x1 = 4 p_x2 = 5 because x12_ref = [1,0,0]
        end

        % System dynamics
        function x_dot = ffct(obj,x_i,u_i)
            G = [cos(x_i(3)), 0; ...
                sin(x_i(3)), 0;
                0, 1];
            x_dot = G * [u_i(1); u_i(2)];
        end

        % Jacobian w.r.t.the states of system dynamics
        function jac_x = dfdx(obj,x_i,u_i)
             jac_x = [0, 0, -u_i(1) * sin(x_i(3));
                 0, 0, u_i(1) * cos(x_i(3));
                 0 0 0];
        end

        function ui = controlInput(obj,x_i,lambda)
            G = [cos(x_i(3)), 0; ...
                sin(x_i(3)), 0; ...
                0, 1];
            u0 = -inv(obj.R) * G' * lambda;
            for i = 1:length(u0)
                if u0(i) < obj.umin(i)
                    u0(i) = obj.umin(i);
                elseif u0(i) > obj.umax(i)
                    u0(i) = obj.umax(i);
                end
            end
            ui = u0;
        end

        % cost function
        function cost = lfct(obj,x_i,u_i)
            cost = 0.5 * (x_i - obj.x_iref)' * obj.Q * (x_i - obj.x_iref) + 0.5 * u_i'*obj.R * u_i;
        end

        % Gradient w.r.t. the states of cost function
        function grad_x = dldx(obj,x_i)
            grad_x = obj.Q * (x_i - obj.x_iref);
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