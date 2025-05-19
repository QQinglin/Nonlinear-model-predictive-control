function  run_sensi(param,agents)
% Performs the iteration of the sensitivity-based Algorithm
num_agents = length(agents);   %{agent_1_sensi,agent_2_sensi,agent_3_sensi};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exchange initial states
for i=1:num_agents
    % exchange required Quantities
    cur_agent = agents{i};   %{agent_1_sensi,agent_2_sensi,agent_3_sensi}; cur_agent = Agent(init,opt_param,agent_model,agent_id);
    num_neighbors = length(cur_agent.neighbors);  % here is {[], []}
    for j=1:num_neighbors
        cur_neighbor = cur_agent.neighbors{j}; %  cur_neighbor = j=1 -->  %  --> obj.neighbors{2} = Neighbor(param, param.A{2, 3}, 1 or 3);
        cur_neighbor_global_ID = cur_neighbor.neighbor_ID;    % neighbor_ID = [1 3] : @isempty get [true, false, true];find get [1 3]
        
        % cur_agent.agent_ID = 2,state = cur_agent.getAgentState() ,
        % agents{1}.setNeighborsStates(cur_agent.getAgentState(),2); 
        agents{cur_neighbor_global_ID}.setNeighborsStates(cur_agent.getAgentState(),cur_agent.agent_ID); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start with q = 1
for q = 1:param.sensi_iter
    
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
        cur_agent = agents{i};      % 2
        num_neighbors = length(cur_agent.neighbors);   % num_neighbors = 2
        for j=1:num_neighbors
            cur_neighbor = cur_agent.neighbors{j}; 
            cur_neighbor_global_ID = cur_neighbor.neighbor_ID;   %  neighbor_ID = [1 3]
            agents{cur_neighbor_global_ID}.setNeighborsStates(cur_agent.getAgentState(),cur_agent.agent_ID);
        end
    end
end
end