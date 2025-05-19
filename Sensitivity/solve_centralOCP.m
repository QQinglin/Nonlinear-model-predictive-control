function  solve_centralOCP(central,u_0,opt_param)
            % solves the OCP via the projected gradient method
            % initial control trajecrtory guess
            u_q = u_0;

            % initialstate, time distribution, step size
            x0 = central.x0;
            t_span = opt_param.t;
            T_end = t_span(end);

            % Convergence
            u_res = zeros(opt_param.grad_iter,1);
            %% Begin of Algorithm
            % first state trajectory
            x_q = ode4(@(t,x)central.dH_dlambda(t,x,u_q,t_span),t_span,x0)';

            for q=1:opt_param.grad_iter

                % backward integration of the adjoint state dynamics
                rev_t_span = flip(t_span);
                x_T = x_q(:,end);
                lambda_T = central.dVdx(T_end,x_T);
                lambda_q_rev = ode4(@(t,lambda)central.dH_dx(t,x_q,u_q,lambda,t_span),rev_t_span,lambda_T);
                lambda_q = flip(lambda_q_rev)';

                % Compute Gradient of OCP w.r.t. input
                grad_u_q = central.dH_du(t_span,x_q,u_q,lambda_q);

                % adapt step size
                alpha_j = opt_param.alpha;

                % projected gradient step
                prev_u_q = u_q;
                u_q = project_box(u_q - alpha_j * grad_u_q,central.umin,central.umax);

                % Forward integration of the state dynamics
                 x_q = ode4(@(t,x)central.dH_dlambda(t,x,u_q,t_span),t_span,x0)';

                % residual
                u_res(q) = compute_residual(u_q, prev_u_q);
            end

            % distribute states, controls adjoint states to agents 
            central.setAgentsStates(x_q, u_q, lambda_q);

end