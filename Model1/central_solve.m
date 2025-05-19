function [u_k,lambda_k, x_k] = central_solve(opt_param, model)
    % Solves the OCP via the fixed point iteration method
    x_0 = model.x0;
    lambda_j = opt_param.lambda0;
    t_span = opt_param.t;  
    T_end = t_span(end);   
    %% Begin of Algorithm
    for j = 1:40
        % forward integration 
         x_j = ode4(@(t, x) model.calcGradHLambda(t, x, lambda_j,t_span), t_span, x_0); % 21 x3

        % backward integration of the adjoint state dynamics
        rev_t_span = flip(t_span);
        x_T = x_j(end, :)';
        lambda_T = model.dVdx(x_T);
        lambda_j_rev = ode4(@(t, lambda) model.calcGradHX(t,lambda,x_j,t_span), rev_t_span, lambda_T);
        lambda_j = flip(lambda_j_rev);  % 21x3
    end
    %% Output
    u_i = zeros(size(x_j,1),size(x_j,2));
    for i = 1:size(x_j, 1)
        u_i(i,:) = model.controlInput(x_j(i,:)', lambda_j(i,:)');
    end
    u_k = u_i;
    x_k = x_j;
    lambda_k = lambda_j;
    %% draw
    figure;
    plot(t_span, x_k(:, 1), '-x', t_span, x_k(:, 2), '-x',t_span, x_k(:, 3), '-x',t_span, x_k(:, 4), '-x',t_span, x_k(:, 5), '-x');
    title('Central Problem Input x_i');
    xlabel('Time');
    ylabel('State Variables');
    legend('x_1', 'x_2','x_3','x_4','x_5');

    figure;
    plot(t_span, lambda_k(:, 1), '-o', t_span, lambda_k(:, 2), '-o',t_span, lambda_k(:, 3), '-o',t_span, lambda_k(:, 4), '-o',t_span, lambda_k(:, 5), '-o');
    title('Central Problem Input \lambda_i');
    xlabel('Time');
    ylabel('Adjoint Variables');
    legend('\lambda_1', '\lambda_2','\lambda_3','\lambda_4','\lambda_5');

    % Optional: Plot control input over time if needed
    figure;
    plot(t_span, u_k(:, 1), '-x', t_span, u_k(:, 2), '-x',t_span, u_k(:, 3), '-x',t_span, u_k(:, 4), '-x',t_span, u_k(:, 5), '-x');
    title('Central Problem Input u_i');
    xlabel('Time');
    ylabel('Control Input');
    legend('u_1', 'u_2','u_3','u_4','u_5');
end