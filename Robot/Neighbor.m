classdef Neighbor < handle
    properties
        x_j;
        lambda_j;
        neighbor_ID{mustBeInteger};
        couplingModel;
        gji_x;
    end

    methods
        function obj = Neighbor(param,model,global_id)

            % global Id for identification
            obj.neighbor_ID = global_id;        

            % neighbor Model 
            obj.couplingModel = model;          

            % initialize neighbor qunatities to receive 
            obj.x_j = zeros(param.N_d,model.nxj);  
            obj.lambda_j = zeros(param.N_d,model.nxj); 
            
            obj.gji_x = zeros(param.N_d,model.nxi);      
         
        end 
    end 
end

