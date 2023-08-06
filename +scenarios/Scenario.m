classdef Scenario < handle
    % SCENARIO 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        name
        x_range             % ROI in x axis
        y_range             % ROI in y axis
        layer_num 
        ISD      
        BS_height  
        Tx_power  
        min_d
    end
    properties(Dependent)
        BSposition
    end
    
    methods
        function obj = Scenario()
            %SCENARIO 构造此类的实例
        end
        
        function BS_pos_list = get.BSposition(obj)
            BS_pos_list = [0 0];
            for layer = 2:obj.layer_num
                vertex = obj.ISD*(layer-1).*exp(1j*[30 90 150 210 270 330 30]./180*pi);
                for n = 1:6
                    d = (vertex(n+1)-vertex(n))/(layer-1);
                    for m = 1:layer-1
                        pos = vertex(n)+d*m;
                        BS_pos_list = [BS_pos_list;[real(pos), imag(pos)]]; %#ok<AGROW>
                    end
                end
            end
            BS_pos_list(abs(BS_pos_list)<10^-7) = 0;
        end
        
        function set.BSposition(obj,layer_num)
            obj.layer_num = layer_num;
            BS_pos_list = [0 0];
            for layer = 2:obj.layer_num
                vertex = obj.ISD*(layer-1).*exp(1j*[30 90 150 210 270 330 30]./180*pi);
                for n = 1:6
                    d = (vertex(n+1)-vertex(n))/(layer-1);
                    for m = 1:layer-1
                        pos = vertex(n)+d*m;
                        BS_pos_list = [BS_pos_list;[real(pos), imag(pos)]]; %#ok<AGROW>
                    end
                end
            end
            BS_pos_list(abs(BS_pos_list)<10^-7) = 0;
            obj.BSposition = BS_pos_list;
        end        
    end
end

