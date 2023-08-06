classdef InF < handle
    % INH_3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        name
        InF_case
        r
        hc
        d_cluster
        d_subsce
        p_subsce
        L
        W
        D
        H
        V
        S
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
        function obj = InF(varargin)
            %INH_3D 构造此类的实例
            %   此处显示详细说明
            obj.name    = 'InF';
            if numel(varargin)>0
                obj.InF_case = varargin{1};
            end
            if strcmp(obj.InF_case,'SL') || strcmp(obj.InF_case,'DH')
                obj.L       = 120;
                obj.W       = 60;
                obj.D       = 20;
            else
                obj.L       = 300;
                obj.W       = 150;
                obj.D       = 50;
            end
            obj.x_range = [-obj.L/2, obj.L/2];
            obj.y_range = [-obj.W/2, obj.W/2];
            if strcmp(obj.InF_case,'SL') || strcmp(obj.InF_case,'SH')
                obj.r  = 0.2;
                obj.hc = 2;
                obj.d_cluster = 10;
                obj.H = 25;
            else
                obj.r  = 0.6;
                obj.hc = 6;
                obj.d_cluster = 2;
                obj.H = 15;
            end
            if strcmp(obj.InF_case,'SL') || strcmp(obj.InF_case,'DL')
                obj.BS_height  = 1.5;  % 8
            else
                obj.BS_height  = 8;  % 8
            end
            if strcmp(obj.InF_case,'DH')
                obj.d_subsce = 1;
                obj.p_subsce = 0.6;
            else
                obj.d_subsce = 0;
                obj.p_subsce = 1;
            end
            obj.V = obj.L*obj.W*obj.H;
            obj.S = 2*(obj.L*obj.W+obj.L*obj.H + obj.W*obj.H);
            obj.min_d      = 1; 
            obj.Tx_power   = 30;    % dB
            obj.BW         = 100e6;
            obj.BS_num     = 18;
            obj.BSposition = [-5/2*obj.D,-obj.D; -3/2*obj.D,-obj.D;-1/2*obj.D,-obj.D;1/2*obj.D,-obj.D;3/2*obj.D,-obj.D;5/2*obj.D,-obj.D;...
                              -5/2*obj.D,0; -3/2*obj.D,0;-1/2*obj.D,0;1/2*obj.D,0;3/2*obj.D,0;5/2*obj.D,0;...
                              -5/2*obj.D,obj.D; -3/2*obj.D,obj.D;-1/2*obj.D,obj.D;1/2*obj.D,obj.D;3/2*obj.D,obj.D;5/2*obj.D,obj.D;...
                ];
        end
        
    end
end

