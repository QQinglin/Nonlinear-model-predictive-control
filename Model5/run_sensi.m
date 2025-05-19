function  run_sensi(param,agents)
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
        % agents{1}.setNeigh(agent2.getAgentS,2) 
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
%% draw
    % draw state figure
    for i = 1:num_agents
        cur_agent = agents{i};
        figure;
        plot(param.t, cur_agent.x_q, '-o');
        hold on;
    
        title(['Agent ', num2str(i), ' xi']);
        xlabel('Time');
        ylabel('State Variables');
        legend('x_1', 'x_2','x_3');
    end
    
    % draw joint variable figure
    for i = 1:num_agents
        cur_agent = agents{i};
        figure;
        plot(param.t, cur_agent.lambda_q, '-o');
        hold on;
    
        title(['Agent ', num2str(i), ' lambda_i']);
        xlabel('Time');
        ylabel('joint Variables');
        legend('lambda_1', 'lambda_2','lambda_3');
    end

    %  % draw input figure
    %     figure;
    % for i = 1:num_agents
    %     cur_agent = agents{i};
    %     plot(param.t, cur_agent.u_q, '-o');
    %     hold on;
    % end
    %     title('Sensi_ui ');
    %     xlabel('Time');
    %     ylabel('State Variables');
    %     legend('u_1', 'u_2','u_3','u_4','u_5');
end