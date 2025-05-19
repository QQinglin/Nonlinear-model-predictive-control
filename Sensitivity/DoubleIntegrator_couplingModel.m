classdef DoubleIntegrator_couplingModel < handle
    % Dynamics of a double integrator

    properties
        nxi{mustBeInteger} = 3;                                   % Number of states
        nui{mustBeInteger} = 3;                                   % Number of controls
        nxj{mustBeInteger} = 3;                                   % Number of states
        modelparams(:,1);
        coupl_fun;
        coupl_jac;
        coupl_cost;
        coupl_cost_jac;
    end

    methods
        function obj = DoubleIntegrator_couplingModel(params,coupl_fun,coupl_jac,coupl_cost,coupl_cost_jac)
            obj.modelparams = params.modelparams;
            obj.coupl_fun = coupl_fun;
            obj.coupl_jac = coupl_jac;
            obj.coupl_cost = coupl_cost;
            obj.coupl_cost_jac = coupl_cost_jac;
        end

        % System dynamics
        function x_dot = fij(obj,~,~)
            %eps = obj.modelparams(1);
            x_dot = obj.coupl_fun;
        end

        % Jacobian w.r.t.the agent states of system dynamics
        function jac_x_i = dfij_dxi(obj,~,~)
            %eps = obj.modelparams(1);
            jac_x_i = obj.coupl_jac;
        end

        % cost function
        function cost = lijfct(~,~)
            cost = obj.coupl_cost;
        end

        % Gradient w.r.t. the states of cost function
        function grad_x_i = dlij_dxi(~,~)
            grad_x_i = obj.coupl_cost_jac;
        end

    end
end