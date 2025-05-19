classdef central_model < handle
    properties
        nx {mustBeInteger} = 3;           % Number of states
        nu {mustBeInteger} = 2;           % Number of controls
        nlambda {mustBeInteger} = 3;      % Number of adjoint variable
        P (3, 3) = ones(3);               % Penalty matrix for the terminal cost
        Q (3, 3) = ones(3);               % Penalty matrix for the state cost
        R (2, 2) = ones(2);               % Penalty matrix for the control input cost
        umin = [];                % Minimum control input
        umax = [];                % Maximum control input
        x0 (3, 1) = ones(3, 1);           % Initial state
        x_iref;
        % x_ijref;


    end

    methods
        function obj = central_model(params)
            obj.P = params.P;
            obj.Q = params.Q;
            obj.R = params.R;
            obj.x0 = params.x0;
            obj.umin = params.umin;
            obj.umax = params.umax;

            obj.x_iref = [3,5,2]';
            % obj.x_ijref = [1,0,0]';
        end

        % System dynamics: compute state derivatives
        function x_dot = ffct(obj, x_i, u_i)
            G = [cos(x_i(3)), 0; ...
                sin(x_i(3)), 0;
                0, 1];
            x_dot = G * [u_i(1); u_i(2)];
        end

        % Jacobian matrix of dynamics
        function jac_x = dfdx(obj,x_i,u_i)
             jac_x = [0, 0, -u_i(1) * sin(x_i(3));
                 0, 0, u_i(1) * cos(x_i(3));
                 0, 0, 0];
        end

        function cost = lfct(obj, x_i, u_i)
            cost_1 = 0.5 * (x_i - obj.x_iref)' * obj.Q * (x_i - obj.x_iref) + 0.5 * u_i'* obj.R * u_i;
            cost_2 = 0 ;% 0.5*(x_j - x_i + obj.x_ijref)' * obj.Q * (x_j - x_i + obj.x_ijref);
            cost = cost_1 + cost_2;
        end

        function grad_x = dldx(obj, x_i)
            grad_x = obj.Q * (x_i - obj.x_iref);
        end

        function cost = Vfct(obj, x_T)
            cost = 0.5 * x_T' * obj.P * x_T;
        end

        function grad_x_V = dVdx(obj, x_T)
            grad_x_V = obj.P * x_T';
        end

        function H = calcHamiltonian(obj, x, u, lambda)
            x_dot = obj.ffct(x, u);
            cost = obj.lfct(x, u);
            H = lambda' * x_dot + cost;
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

        function grad_h_lambda = calcGradHLambda(obj,t, x, lambda,t_span)
             % Compute the gradient of Hamiltonian w.r.t. lambda
             lambda_t = interp1(t_span', lambda, t)'; % t_span 1x21, lambda 21 x 3, x: 3x1, lambda_t 1x3
             u = obj.controlInput(x,lambda_t);

             grad_h_lambda = obj.ffct(x, u);
        end

        function grad_h_x = calcGradHX(obj,t,x_q,lambda,t_span)
            x = interp1(t_span',x_q,t)'; % 3x1
            u = obj.controlInput(x,lambda);

            grad_h_x = -obj.dfdx(x, u)'*lambda - obj.dldx(x);
         end
    end
end
