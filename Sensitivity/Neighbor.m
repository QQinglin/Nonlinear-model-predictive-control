classdef Neighbor < handle
    % Takes care of all the local steps

    properties
        x_j;
        lambda_j;
        neighbor_ID{mustBeInteger};
        couplingModel;
        % gradients 
        gji_x;
    end

    methods
        function obj = Neighbor(param,model,global_id)
            % global Id for identification
            obj.neighbor_ID = global_id;

            % neighbor Model 
            obj.couplingModel = model;

            % initialize neighbor qunatities 
            obj.x_j = zeros(model.nxj,param.N_d); 
            obj.lambda_j = zeros(model.nxj,param.N_d); 

           obj.gji_x = zeros(model.nxi,param.N_d); 
        end 
    end 
end

