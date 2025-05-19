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

% Start with q = 1
for q = 1:param.sens_iter_k

    % Solution of local OCP with sensitivities
    for i=1:num_agents
        % Compute Sensitivities
        agents{i}.compute_sensitivities();
        % Solve local OCP
        agents{i}.solve_localOCP();

        % Store the delta values in the 3D arrays
        delta_kx(q,i) = max(abs(x_k(:,i) - agents{i}.x_q));
        delta_klambda(q,i) = max(abs(lambda_k(:,i) - agents{i}.lambda_q));
        delta_ku(q,i) = max(abs(u_k(:,i) - agents{i}.u_q));
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
for i = 1:num_agents
        cur_agent = agents{i};
        for j = 1:param.sens_iter_k - 1
            cur_agent.collect_result = [cur_agent.collect_result, max(cur_agent.first_x(:,j+1) - cur_agent.last_x(:,j))];
        end

        matrix_1 = cur_agent.collect_result;
        
        figure;  
        plot(1:param.sens_iter_k-1, matrix_1, '-o');
        hold on;

        title(['Evaluation ', num2str(i), ' xi']);
        xlabel('Iteration k');
        ylabel('$$\|\delta x_i^k\|_{L_\infty}$$', 'Interpreter', 'latex');
        legend('||\delta x_1^k_{L_\infty}||');
        hold off;  
 end

%% draw
    data_eva = [linspace(0, param.sens_iter_k, 30)]'; % *****************
    data_eva = [data_eva, delta_kx]; % **************
    writematrix(data_eva, 'C:\Users\24966\Desktop\project\3eva.csv'); % *********************

    figure;
    for i = 1:size(delta_kx, 2)
        plot(linspace(0, param.sens_iter_k, 30), delta_kx(:, i), 'x-'); % Use 'x-' to draw lines connecting the points
        hold on;
    end
    title('Evaluation');
    xlabel('Sensitivity Iteration');
    ylabel('X*-X');
    legend('x_1', 'x_2', 'x_3');
    hold off;

    % draw state figure
    data_xq = [param.t']; % *****************

    figure;
    for i = 1:num_agents
        cur_agent = agents{i};
        data_xq = [data_xq, cur_agent.x_q]; % **************
        plot(param.t, cur_agent.x_q, '-o');
        hold on;
    end

        writematrix(data_xq, 'C:\Users\24966\Desktop\project\3Sen_x_q.csv'); % *********************
        title('Sensi x_i ');
        xlabel('Time');
        ylabel('State Variables');
        legend('x_1', 'x_2','x_3');

    % draw joint variable figure
    data_lambdaq = [param.t']; % *****************

    figure;
    for i = 1:num_agents
        cur_agent = agents{i};
        data_lambdaq = [data_lambdaq, cur_agent.lambda_q]; % **************
        plot(param.t, cur_agent.lambda_q, '-o');
        hold on;
    end
        writematrix(data_lambdaq, 'C:\Users\24966\Desktop\project\3Sen_lambda_q.csv'); % *********************
        title('Sensi \lambda_i ');
        xlabel('Time');
        ylabel('joint Variables');
        legend('\lambda_1', '\lambda_2','\lambda_3');

     % draw input figure
        data_uq = [param.t']; % *****************

        figure;
    for i = 1:num_agents
        cur_agent = agents{i};
        data_uq = [data_uq, cur_agent.u_q]; % **************
        plot(param.t, cur_agent.u_q, '-o');
        hold on;
    end
        writematrix(data_uq, 'C:\Users\24966\Desktop\project\3Sen_u_q.csv'); % *********************
        title('Sensi u_i ');
        xlabel('Time');
        ylabel('State Variables');
        legend('u_1', 'u_2','u_3');    
end