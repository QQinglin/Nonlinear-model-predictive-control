function  run_sensi(param,agents,u_k,lambda_k,x_k)
% Performs the iteration of the sensitivity-based Algorithm
num_agents = length(agents);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exchange initial states
for i=1:num_agents
    % exchange required Quantities
    cur_agent = agents{i};        % agents = {agent_1_sensi,agent_2_sensi,agent_3_sensi}      
    num_neighbors = length(cur_agent.neighbors);  % {[obj.neigh_1],[obj.neigh_3]}
    for j=1:num_neighbors
        cur_neighbor = cur_agent.neighbors{j};
        cur_neighbor_global_ID = cur_neighbor.neighbor_ID;
        % agents{1}.setNeighS(agent2.getAgentS,2) 
        agents{cur_neighbor_global_ID}.setNeighborsStates(cur_agent.getAgentState(),cur_agent.agent_ID);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 3;
tol_res = 1e-2;

% Start with q = 1
for k = 1:param.sens_iter_k

    % Solution of local OCP with sensitivities
    for i=1:num_agents
        % Compute Sensitivities
        agents{i}.compute_sensitivities();
        % Solve local OCP
        agents{i}.solve_localOCP();

        %%
        % Store the delta values in the 3D arrays
        delta_kx(k,i) = max(abs(x_k(:,i) - agents{i}.x_q));
        delta_klambda(k,i) = max(abs(lambda_k(:,i) - agents{i}.lambda_q));
        delta_ku(k,i) = max(abs(u_k(:,i) - agents{i}.u_q));
    %% Anderson acceleration part
        % agents{i}.sen_k = [agents{i}.sen_k,abs(agents{i}.x_q - agents{i}.last_re(:,end))];
        agents{i}.sen_k = [agents{i}.sen_k,max(abs(agents{i}.x_q - agents{i}.last_re(:,end)))];
        % if max(abs(agents{i}.x_q - agents{i}.last_re(:,end))) < tol_res
        %     agents{i}.sen_k = [agents{i}.sen_k,k]; % record all convergence iteration
        % end

        agents{i}.last_re = [agents{i}.last_re, agents{i}.x_q]; % record the last iteration result 
        %%
        if k == 1 
            agents{i}.x_bar = [agents{i}.x_bar, agents{i}.x_q]; % x_bar (1)
            agents{i}.x = [agents{i}.x,agents{i}.x_bar(:,end)]; % x(2) = x1, x(1) = x0
        end

        agents{i}.conver_k = [agents{i}.conver_k,max(abs(agents{i}.x(:,end) - agents{i}.x(:,size(agents{i}.x,2) - 1)))];
        % if max(abs(agents{i}.x(:,end) - agents{i}.x(:,size(agents{i}.x,2) - 1))) < tol_res
        %     agents{i}.conver_k = [agents{i}.conver_k,k]; % record all convergence iteration after using Anderson
        % end

        if k >= 2
            m_k = min(k-1, m); % in the paper k = 1 at beginning, here k = 2 at beginning, so k-1
            agents{i}.x_bar = [agents{i}.x_bar, agents{i}.x_q]; % x_bar (2)

            fun = @(alpha) 0; % creating a anonymous function 
            for j = 1 : m_k+1 % m_k = 1 , m = 3
                fun = @(alpha) fun(alpha) + ...
                alpha(m_k - j + 2)^2 * sum((agents{i}.x_bar(:,k - j + 1) - agents{i}.x(:,k - j + 1)).^2);
            end

            alpha0 = ones(m_k + 1, 1) / (m_k + 1); % need initial points for using fmincon command
            Aeq = ones(1, m_k+1); % the sum of alpha = 1
            beq = 1;

            alpha = fmincon(fun, alpha0, [], [], Aeq, beq);

            % X(k+1) = alpha * x_bar
            agents{i}.x = [agents{i}.x, sum(agents{i}.x_bar(:,size(agents{i}.x_bar,2) - m_k : end) .* alpha', 2)];
        end  
    
    end

    % Communication Step
    for i=1:num_agents
        % exchange required Quantities
        cur_agent = agents{i};
        num_neighbors = length(cur_agent.neighbors);
        for j=1:num_neighbors
            cur_neighbor = cur_agent.neighbors{j};
            cur_neighbor_global_ID = cur_neighbor.neighbor_ID;
            agents{cur_neighbor_global_ID}.setNeighborsStates(cur_agent.getAgentState(),cur_agent.agent_ID);
        end
    end
end

%%
 anderson_sen_x = [1:param.sens_iter_k]';
 no_anderson_sen_x = [1:param.sens_iter_k]';
 for i = 1:num_agents
        cur_agent = agents{i};
        index_sen = find(cur_agent.sen_k < 1e-2, 1);
        index_cen = find(cur_agent.conver_k < 1e-2, 1);

        fprintf ('no Anderson conver_k = %f\n', index_sen);
        fprintf ('with Anderson conver_k = %f\n', index_cen);

        anderson_sen_x = [anderson_sen_x,cur_agent.sen_k'];
        no_anderson_sen_x = [no_anderson_sen_x,cur_agent.conver_k'];

        yline(1e-2, '--r', '1e-2 Threshold', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');

        figure;  
        
        plot(1:param.sens_iter_k, cur_agent.sen_k, '-o');
        hold on;
        plot(1:param.sens_iter_k, cur_agent.conver_k, '-x'); % '-x' 表示使用不同的标记
        title(['Evaluation ', num2str(i), ' xi']);
        hold off;

         writematrix(anderson_sen_x, 'C:\Users\24966\Desktop\project\anderson_xq.csv');
         writematrix(no_anderson_sen_x, 'C:\Users\24966\Desktop\project\noanderson_xq.csv');
 end


%%
% for i = 1:num_agents
%         cur_agent = agents{i};
%         for j = 1:param.sens_iter_k - 1
%             cur_agent.collect_result = [cur_agent.collect_result, max(cur_agent.first_x(:,j+1) - cur_agent.last_x(:,j))];
%         end
% 
%         matrix_1 = cur_agent.collect_result;
% 
%         figure;  
%         plot(1:param.sens_iter_k-1, matrix_1, '-o');
%         hold on;
% 
%         title(['Evaluation ', num2str(i), ' xi']);
%         xlabel('Iteration k');
%         ylabel('$$\|\delta x_i^k\|_{L_\infty}$$', 'Interpreter', 'latex');
%         legend('||\delta x_1^k_{L_\infty}||');
%         hold off;  
%  end

    % for i = 1:num_agents
    %     fprintf ('no Anderson conver_k = %f\n', agents{i}.conver_k(1));
    % 
    %     fprintf('sen_k = %f\n', agents{i}.sen_k(1));
    % 
    % end
%%
% data_eva = [linspace(0, param.sens_iter_k, 30)]'; % *****************
%     data_eva = [data_eva, delta_kx]; % **************
%     writematrix(data_eva, 'C:\Users\24966\Desktop\project\1eva.csv'); % *********************
% 
%     figure;
%     for i = 1:size(delta_kx, 2)
%         plot(linspace(0, param.sens_iter_k, 30), delta_kx(:, i), 'x-'); % Use 'x-' to draw lines connecting the points
%         hold on;
%     end
%     title('Evaluation');
%     xlabel('Sensitivity Iteration');
%     ylabel('X*-X');
%     legend('x_1', 'x_2', 'x_3','x_4','x_5');
%     hold off;
% % Draw state figure
% data_xq = [param.t'];
% figure;
% for i = 1:num_agents
%     cur_agent = agents{i};
%     data_xq = [data_xq, cur_agent.x_q];
%     plot(param.t, cur_agent.x_q, '-o');
%     hold on;
% end
% Save data as CSV
% writematrix(data_xq, 'C:\Users\24966\Desktop\project\xq_data.csv');

% title('Based-Sensitivity x_i');
% xlabel('Time');
% ylabel('State Variables');
% legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5');
% 
% % Draw joint variable figure
% figure;
% 
% data_lambda_q = [param.t'];
% for i = 1:num_agents
%     cur_agent = agents{i};
%     data_lambda_q = [data_lambda_q, cur_agent.lambda_q];
% 
%     plot(param.t, cur_agent.lambda_q, '-o');
%     hold on;
% end
% % Save data as CSV
% writematrix(data_lambda_q, 'C:\Users\24966\Desktop\project\lambdaq_data.csv');
% 
% title('Based-Sensitivity \lambda_i');
% xlabel('Time');
% ylabel('adjoint Variables');
% legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5');
% 
% % Draw input figure
% figure;
% data_u_q = [param.t'];
% for i = 1:num_agents
%     cur_agent = agents{i};
%     data_u_q = [data_u_q, cur_agent.u_q];
%     plot(param.t, cur_agent.u_q, '-o');
%     hold on;
% end
% % Save data as CSV
% writematrix(data_u_q, 'C:\Users\24966\Desktop\project\uq_data.csv');
% 
% title('Based-Sensitivity u_i');
% xlabel('Time');
% ylabel('State Variables');
% legend('u_1', 'u_2', 'u_3', 'u_4', 'u_5');
