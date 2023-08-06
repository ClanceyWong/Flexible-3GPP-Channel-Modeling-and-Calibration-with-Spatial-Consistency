classdef RMa_3D < scenarios.Scenario
    % 3D_RMA 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        W                   % building width
        h                   % building height
        BW                  % bandwidth
    end
    properties(Dependent)
        R   % radius
    end
    
    methods
        function obj = RMa_3D(varargin)
            %3D_RMA 构造此类的实例
            if ~isempty(varargin)
                Layer_num  = varargin{1};
            else
                Layer_num  = 3;
            end
            obj.name       = '3D-RMa';
            obj.W          = 20;
            obj.h          = 5;
            obj.ISD        = 1000;            % 500 m (option: 200m) for 3D-UMa; 200 m for 3D-UMi.
            obj.BS_height  = 35;             % 25 m for 3D-UMa; 10 m for 3D-UMi.
            obj.Tx_power   = 46;             % dBm
            obj.min_d      = 35;             % the min 2D distance between BS and UE, 35 m for 3D-UMa; 10 m for 3D-UMi.
            obj.layer_num  = Layer_num;      % the number of layer of cells.
            obj.x_range    = [-ceil(obj.R*(1.5*obj.layer_num-0.5)), ceil(obj.R*(1.5*obj.layer_num-0.5))];
            obj.y_range    = [-ceil(obj.ISD*(obj.layer_num-0.5)), ceil(obj.ISD*(obj.layer_num-0.5))];
            obj.BW         = 10e6;
        end
        
        function value = get.R(obj)
            value = obj.ISD/sqrt(3);
        end
    end
end
