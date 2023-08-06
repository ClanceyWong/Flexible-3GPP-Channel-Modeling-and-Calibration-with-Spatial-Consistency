classdef BS_Sector < handle
    %BS_SECTOR 
    
    properties
        ID           % the ID of sector [attached_BS.ID, local ID, global ID]
        frequency    % the RF in Hz
        BW           % the bandwidth
        boresight    % the direction of the sector
        antenna      % the configurational antenna
        attached_BS  % the attached BaseStation of current sector
        attached_UE  % the ID of attached UE
        
        PHI_n_m
        link_params
        
        lineColor    % used to plot the links
    end
    properties(Dependent)
        num_UE       % the number of attached UE
    end
    
    methods
        function obj = BS_Sector(varargin)
            %BS_SECTOR 
            if length(varargin) >= 1
                obj.attached_BS = varargin{1};
            end
            if length(varargin) == 2
                obj.frequency = varargin{2};
            end
            if isempty(obj.lineColor)
                obj.lineColor = 0.8*rand(1,3)+0.1;
            end
        end
        
        function addUE(obj,newUE)
            tmp = [obj.attached_UE, newUE];
            obj.attached_UE = tmp;
        end
%         function UElist = getUE(obj,total_list)
%             UElist = total_list(obj.attached_UE);
%         end
        
        
        function init_link_params(obj,link)
            if isempty(obj.link_params)
                indx = 1;
            else
                indx = length(obj.link_params)+1;
            end
            obj.link_params{indx}.UEID              = link.UE.ID;
            obj.link_params{indx}.BS_pos_wrap       = link.BS_pos_wrap;
            obj.link_params{indx}.O2I               = link.O2I;
            obj.link_params{indx}.bLOS              = link.bLOS;
            obj.link_params{indx}.K                 = link.K;        % rice factor [1x1]
            obj.link_params{indx}.tau_n             = link.tau_n;    % cluster delay [1xN]
            obj.link_params{indx}.tau_n_LOS         = link.tau_n_LOS;% cluster delay [1xN]
            obj.link_params{indx}.Pn                = link.Pn;       % cluster power [1xN]
            obj.link_params{indx}.strong_cluster_id = link.strong_cluster_id;
            obj.link_params{indx}.map_delay         = link.map_delay;
            obj.link_params{indx}.theta_n_m_ZOA     = link.theta_n_m_ZOA;
            obj.link_params{indx}.theta_n_m_ZOD     = link.theta_n_m_ZOD;
            obj.link_params{indx}.phi_n_m_AOA       = link.phi_n_m_AOA;
            obj.link_params{indx}.phi_n_m_AOD       = link.phi_n_m_AOD;
            obj.link_params{indx}.theta_LOS_ZOA     = link.theta_LOS_ZOA;
            obj.link_params{indx}.theta_LOS_ZOD     = link.theta_LOS_ZOD;
            obj.link_params{indx}.phi_LOS_AOA       = link.phi_LOS_AOA;
            obj.link_params{indx}.phi_LOS_AOD       = link.phi_LOS_AOD;
            obj.link_params{indx}.XPR_n_m           = link.XPR_n_m;
            obj.link_params{indx}.PL                = link.PL;
            obj.link_params{indx}.SF                = link.SF;
        end
        
        function plot_links(obj)
            pos = zeros(6,obj.num_UE);
            pos(1:2:5,:) = [obj.attached_UE(:).pos3D];
            pos(2:2:6,:) = [reshape([obj.attached_UE(:).attachBS_pos],2,[]); obj.attached_BS.h_BS*ones(1,obj.num_UE)];
            figure(1); hold on;
            plot3(pos(1:2,:),pos(3:4,:),pos(5:6,:),'color',obj.lineColor); hold on;
        end
        
        %% get and set function
        function out = get.num_UE(obj)
            out = length(obj.attached_UE);
        end
    end
end

