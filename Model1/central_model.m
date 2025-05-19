classdef central_model < handle
    properties
        nx {mustBeInteger} = 5;           % Number of states
        nu {mustBeInteger} = 5;           % Number of controls
        nlambda {mustBeInteger} = 5;      % Number of adjoint variable
        P ;               % Penalty matrix for the terminal cost
        Q ;               % Penalty matrix for the state cost
        R ;               % Penalty matrix for the control input cost
        umin ;               % Minimum control input
        umax ;                % Maximum control input
        x0 ;          % Initial state
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
            x_1 = x(2) + x(1) * u(1);
            x_2 = x(1) + x(3) + x(2) * u(2) + x(5) - x(2);
            x_3 = x(4) + x(2) +  x(3) + x(3) * u(3);
            x_4 = x(3) + 2 * x(4) + x(5) +u(4) * x(4);
            x_5 = x(2) + x(4) + 3 * x(5) + u(5) * x(5);
            % x_dot = [x(2) + x(1) * u(1);
            %     x(1) + x(3) + x(2) * u(2) + x(5) - x(2);
            %     x(4) + x(2) +  x(3) + x(3) * u(3);
            %     x(3) + 2 * x(4) + x(5) +u(4) * x(4);
            %     x(2) + x(4) + 3 * x(5) + u(5) * x(5)];
            x_dot = [x_1;x_2;x_3;x_4;x_5];

        end

        % Jacobian matrix of dynamics
        function f_dot = dfdx(obj, x, u)
            f_dot = [u(1), 1, 0, 0, 0;
                     1, u(2)-1, 1, 0, 1;
                     0, 1, 1 + u(3),1 ,0;
                     0, 0, 1, 2+u(4), 1;
                     0, 1, 0, 1, 3+u(5)];
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
            % g_x = [x(1),0,0,0,0;
            %     0,x(2),0,0,0;
            %     0,0,x(3),0,0;
            %     0,0,0,x(4),0;
            %     0,0,0,0,x(5)];
             u0 = -inv(obj.R) * lambda .* x; % R 3x3
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
             lambda_t = interp1(t_span, lambda, t)'; 
             u = obj.controlInput(x,lambda_t);

            grad_h_lambda = obj.ffct(x, u);
        end

        function grad_h_x = calcGradHX(obj, t,lambda,x,t_span)
            x = interp1(t_span',x,t)'; % 3x1
            u = obj.controlInput(x,lambda);

            grad_h_x = -obj.dfdx(x, u)'*lambda - obj.dldx(x);
         end
    end
end
