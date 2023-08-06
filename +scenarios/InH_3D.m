classdef InH_3D < handle
    % INH_3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        name
        InH_case
        BS_num
        BSposition
        BS_height  
        Tx_power  
        BW                  % bandwidth
        min_d
        x_range             % ROI in x axis
        y_range             % ROI in y axis
    end
    
    methods
        function obj = InH_3D(varargin)
            %INH_3D 构造此类的实例
            %   此处显示详细说明
            obj.name    = '3D-InH';
            obj.x_range = [-60, 60];
            obj.y_range = [-25, 25];
            obj.min_d       = 3; 
            obj.BS_height   = 3;
            obj.Tx_power    = 24;    % dB
            obj.BW          = 10e6;
            if ~isempty(varargin)
                obj.BS_num   = varargin{1};
                if numel(varargin)>1
                    obj.InH_case = varargin{2};
                end
            else
                obj.BS_num = 2;
            end
            if obj.BS_num == 2
                obj.BSposition = [-30 0;30 0];
            elseif obj.BS_num == 3
                obj.BSposition = [0 0; -40 0;40 0];
            elseif obj.BS_num == 6 && strcmp(obj.InH_case,'A')
                obj.BSposition = [-10 0; 10 0; -30 0; 30 0; -50 0;50 0];
            elseif obj.BS_num == 6 && strcmp(obj.InH_case,'B')
                obj.BSposition = [0 10; -40 10;40 10; 0 -10; -40 -10;40 -10];
            elseif obj.BS_num == 12 
                obj.BSposition = [-10 10; 10 10; -30 10; 30 10; -50 10;50 10; -10 -10; 10 -10; -30 -10; 30 -10; -50 -10;50 -10];
            end
        end
        
    end
end

