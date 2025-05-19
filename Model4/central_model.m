classdef central_model < handle
    properties
        nx {mustBeInteger} = 3;           % Number of states
        nu {mustBeInteger} = 3;           % Number of controls
        nlambda {mustBeInteger} = 3;      % Number of adjoint variable
        P (3, 3) = ones(3);               % Penalty matrix for the terminal cost
        Q (3, 3) = ones(3);               % Penalty matrix for the state cost
        R (3, 3) = ones(3);               % Penalty matrix for the control input cost
        umin (1, 1) = -inf;               % Minimum control input
        umax (1, 1) = inf;                % Maximum control input
        x0 (3, 1) = ones(3, 1);          % Initial state
        modelparams(:,1);
    end

    methods
        function obj = central_model(params)
            obj.P = params.P;
            obj.Q = params.Q;
            obj.R = params.R;
            obj.x0 = params.x0;
            obj.umin = params.umin;
            obj.umax = params.umax;
        end

        % System dynamics: compute state derivatives
        function x_dot = ffct(obj, x, u)
            x_dot = [x(2) + x(1) * u(1);
                x(1) + x(3) + x(2) * u(2);
                3 * x(3) + x(3) * u(3)];
        end

        % Jacobian matrix of dynamics
        function f_dot = dfdx(obj, x, u)
            f_dot = [u(1), 1, 0;
                     1, u(2), 1;
                     0, 0, 3 + u(3)];
        end

        function cost = lfct(obj, x, u)
            cost = 0.5 * x' * obj.Q * x + 0.5 * u' * obj.R * u;
        end

        function grad_x = dldx(obj, x)
            grad_x = obj.Q * x;
        end

        function cost = Vfct(obj, x_T)
            cost = 0.5 * x_T' * obj.P * x_T;
        end

        function grad_x_V = dVdx(obj, x_T)
            grad_x_V = obj.P * x_T;
        end

        function H = calcHamiltonian(obj, x, u, lambda)
            x_dot = obj.ffct(x, u);
            cost = obj.lfct(x, u);
            H = lambda' * x_dot + cost;
        end

        function ui = controlInput(obj, x, lambda)
             u0 = -inv(obj.R) * lambda .* x; % R 3x3
            if u0 < obj.umin
                ui = obj.umin;
            elseif u0 > obj.umax
                ui = obj.umax;
            else
                ui = u0;
            end
        end

        function grad_h_lambda = calcGradHLambda(obj,t, x, lambda,t_span)
             % Compute the gradient of Hamiltonian w.r.t. lambda
             lambda_t = interp1(t_span, lambda, t)'; % t_span 1x21, lambda 21 x 3, x: 3x1, lambda_t 1x3
             u = obj.controlInput(x,lambda_t);

            grad_h_lambda = obj.ffct(x, u);
        end

        function grad_h_x = calcGradHX(obj, t,lambda,x,t_span)
            x = interp1(t_span,x,t)'; % 3x1
            u = obj.controlInput(x,lambda);

            grad_h_x = -obj.dfdx(x, u)'*lambda - obj.dldx(x);
         end
    end
end
