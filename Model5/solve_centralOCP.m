function  solve_centralOCP(model_param, model)
            % solves the OCP via the fixed point iteration method
            x0 = model.x0; % x0 = [1;1;1]
            lambda_q = model_param.lambda0; % 21 x 3 ones
            t_span = model_param.t;
            T_end = t_span(end);
            %% Begin of Algorithm
            for q=1:model_param.fixed_iter_j
                % forward integration of the state dynamics
                x_q = ode4(@(t,x)model.calcGradHLambda(t,x,lambda_q,t_span),t_span,x0);

                % backward integration of the adjoint state dynamics
                rev_t_span = flip(t_span);
                x_T = x_q(end,:);             % 1x3
                lambda_T = model.dVdx(x_T); % 3x1
                lambda_q_rev = ode4(@(t,lambda)model.calcGradHX(t,x_q,lambda,t_span),rev_t_span,lambda_T);
                lambda_q = flip(lambda_q_rev);
            end

            % distribute states, controls adjoint states to agents 
            % update the output
            u_q = [];
            for j = 1:size(x_q,1)
                u_q = [u_q, model.controlInput(x_q(j,:)',lambda_q(j,:)')]; 
            end

    figure;
    plot(t_span, x_q(:,1), '-o', t_span, x_q(:,2), '-x',t_span, x_q(:,3), '-x');
    title('State Variables over Time');
    xlabel('Time');
    ylabel('State Variables');
    legend('x_1', 'x_2','x_3');

    % figure;
    % plot(t_span, lambda_q(:,3), '-o', t_span, lambda_q(:,2), '-x',t_span, lambda_q(:,3), '-x');
    % title('Adjoint Variables over Time');
    % xlabel('Time');
    % ylabel('Adjoint Variables');
    % legend('lambda_1', 'lambda_2','lambda_3');

    % Optional: Plot control input over time if needed
    % figure
    % plot(t_span, u_q(:,1), '-o', t_span, u_q(:,2), '-x',t_span, u_q(:,3), '-x');
    % title('Input over Time');
    % xlabel('Time');
    % ylabel('Input Variables');
    % legend('u_1', 'u_2','u_3');

    %%
    % opt_param.grad_iter = 400;  
% agent_1_central = agentModel_1; %Agent(init,opt_param,agent_model,1);
% agent_2_central = agentModel_2; %Agent(init,opt_param,agent_model,2);
% agent_3_central = agentModel_3; %Agent(init,opt_param,agent_model,3);
% agents_central  = {agent_1_central,agent_2_central,agent_3_central};
% central = CentralProblem(agents_central);
% solve_centralOCP(central,opt_param);
end