classdef DoubleIntegrator_agentModel < handle
    % Dynamics of a double integrator

    properties
        nx_i {mustBeInteger} = 1;                                 % Number of states
        nu_i {mustBeInteger} = 1;                                 % Number of controls
        P(1,1) double = ones(1);                                % Adjusted size for P
        Q(1,1) double = ones(1);                                % Adjusted size for Q
        R(1,1) double = ones(1);                                % Adjusted size for R
        umin(1,1) double = -inf(1,1);                           % Adjusted size for umin
        umax(1,1) double = inf(1,1);                            % Adjusted size for umax
        x0(1,1) double = ones(1,1);                             % Initial state
        modelparams(:,1) double; 
        dynamic_fun;
        dynamic_jac;

        nx_j{mustBeInteger} = 3;
    end

    methods
        function obj = DoubleIntegrator_agentModel(params,dynamic_fun,dynamic_jac)

            obj.P = params.P;
            obj.Q = params.Q;
            obj.R = params.R;
            obj.x0 = params.x0;
            obj.umin = params.umin;
            obj.umax = params.umax;
            obj.dynamic_fun = dynamic_fun;
            obj.dynamic_jac = dynamic_jac;
            
        end

        % System dynamics
        function x_dot = ffct(x_i,u_i)
            x_dot = obj.dynamic_func;
        end

        % Jacobian w.r.t.the states of system dynamics
        function jac_x = dfdx(x_i,u_i)
            jac_x = obj.dynamic_jac;   %0, 2+u_i(2), 0;0, 0, 3+u_i(3)];
        end


        % cost function
        function cost = lfct(obj,x_i,u_i)
            cost = 0.5 * x_i'*obj.Q*x_i + 0.5 * u_i'*obj.R*u_i;
        end

        % Gradient w.r.t. the states of cost function
        function grad_x = dldx(obj,~,x_i,~)
            grad_x = obj.Q*x_i;
        end

        % Gradient w.r.t. the inputs of cost function
        function grad_u = dldu(obj,~,~,u_i)
            grad_u = obj.R*u_i;
        end

        % terminal cost
        function cost = Vfct(obj,~,xi_T)
            cost = 0.5 * xi_T'*obj.P*xi_T;
        end

        % Gradient w.r.t. to the states of terminal cost
        function grad_x = dVdx(obj,~,xi_T)
            grad_x = obj.P*xi_T;
        end
    end
end