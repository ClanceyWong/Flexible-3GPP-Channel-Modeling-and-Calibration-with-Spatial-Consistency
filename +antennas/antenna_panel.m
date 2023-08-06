classdef antenna_panel < handle
    % ANTENNA_ARRAY
    
    properties
        ID                    % ID
        M                     % M rows
        N                     % N columns
        Kv                    % Table 7.1-1: (Complex weight for antenna element m in elevation)
        Kh                    
        d_H                   % Horizontal antenna element spacing
        d_V                   % Vertical antenna element spacing 
        P                     % 1: vertically polarized; 2: cross-polarized
        X_pol                 % represent the angle of cross-polarized antennas
        ele_downtilt          % electric downtilt
        ele_panning
        
        attachedArray         % the calss of {antenna_array}
        element_list          % the vector of all {antenna_element} class
    end
    properties(Dependent)
        alpha                 % bearing angle for LCS to GCS
        beta                  % downtilt angle for LCS to GCS
        gamma                 % slant angle for LCS to GCS
        attachedType          % 'BS' or 'UE'
        num_element           % the number of antenna element
        pos_LCS               % the position of current antenna panel in LCS
        pos_element_LCS       % the position of antenna element in LCS
        w_m                   % Complex weight for antenna element m in elevation
        w_n                   % Complex weight for antenna element m in elevation
        num_port
    end
    
    methods
        function obj = antenna_panel(varargin)
            % (M=1, N=1, K=1, d_H=0.5, d_V=0.5, P=2, X_pol=[-45, 45], ele_downtilt=90, alpha=0, beta=0, gamma=0,attachedType='BS')
            if isempty(varargin)
                obj.M            = 2; 
                obj.N            = 2; 
                obj.Kv           = 1; 
                obj.Kh           = 1; 
                obj.d_H          = 0.5; 
                obj.d_V          = 0.5; 
                obj.P            = 1; %        obj.P            = 2; 
                obj.X_pol        = 0; %        obj.X_pol        = [-45,45]; 
                obj.ele_downtilt = 102; 
                obj.ele_panning  = 0;
            else
                obj.M            = varargin{1}.M; 
                obj.N            = varargin{1}.N; 
                obj.Kv           = varargin{1}.Kv; 
                obj.Kh           = varargin{1}.Kh; 
                obj.d_H          = varargin{1}.d_H; 
                obj.d_V          = varargin{1}.d_V; 
                obj.P            = varargin{1}.P; 
                obj.X_pol        = varargin{1}.X_pol; 
                obj.ele_downtilt = varargin{1}.ele_downtilt;
                obj.ele_panning  = varargin{1}.ele_panning;
            end
            if obj.P == 1
                obj.X_pol        = 0;
                obj.element_list = antennas.antenna_element(obj.X_pol(1));
                obj.element_list.attachedPanel = obj;
                
                if ~isempty(varargin) && isfield(varargin{1},'pol_model')
                    obj.element_list.pol_model = varargin{1}.pol_model;
                end
            elseif obj.P == 2
                obj.element_list = [antennas.antenna_element(obj.X_pol(1)), antennas.antenna_element(obj.X_pol(2));];
                obj.element_list(1).attachedPanel = obj;
                obj.element_list(2).attachedPanel = obj;
                if ~isempty(varargin) && isfield(varargin{1},'pol_model')
                    obj.element_list(1).pol_model = varargin{1}.pol_model;
                    obj.element_list(2).pol_model = varargin{1}.pol_model;
                end
            end
            
        end
        
        function gain = panel_gain(obj,link)% need to be modified.
            lambda      = 3e8/link.sector.frequency;
            pos_UE      = [link.UE.pos, link.UE.h_UT]';

            Nc        = obj.N;
            pos_BS      = [link.BS_pos_wrap, link.sector.h_BS]';
            r_tx_LOS    = [sind(link.theta_LOS_ZOD)*cosd(link.phi_LOS_AOD);  % the spherical unit vector
                           sind(link.theta_LOS_ZOD)*sind(link.phi_LOS_AOD); 
                           cosd(link.theta_LOS_ZOD)];
            gain_LOS    = zeros(1,Nc);
            PHI_LOS     = 2*pi/lambda*sqrt(sum((pos_UE- pos_BS).^2));
            F_rx_theta_phi = link.UE.antenna.panel(1).element_list(1).field_pattern(link.phi_LOS_AOA, link.theta_LOS_ZOA);
            for s =1:Nc
                k      = 1:obj.Kv;
%                     ID = reshape([link.sector.antennaArray.element_list(k).ID],2,[])';
%                     d_tx_s = [zeros(length(k),1), (ID(:,2)-1)*link.sector.antennaArray.d_H, (ID(:,1)-1)*link.sector.antennaArray.d_V]'*lambda + pos_BS;
                d_tx_s = obj.attachedArray.R*(obj.pos_LCS + permute(obj.pos_element_LCS(k,s,:),[3,1,2]))*lambda;
                F_tx_theta_phi_tmp = obj.element_list(1).field_pattern(link.phi_LOS_AOD, link.theta_LOS_ZOD);
                F_tx_theta_phi     = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*obj.w_m)*F_tx_theta_phi_tmp;
                gain_LOS(1,s)      = sum(abs(F_rx_theta_phi.'*[exp(-1j*PHI_LOS),0;0,-exp(-1j*PHI_LOS)]*F_tx_theta_phi).^2);
            end
            gain = 10*log10(gain_LOS);
        end
        
        function set_ele_downtilt(obj,new_ele_downtilt)
            obj.ele_downtilt = new_ele_downtilt;
        end
        function [indr, indc] = get_port_pos(obj,port_indx)
            if port_indx<0 && port_indx>=obj.num_port
                error(['the port ',num2str(port_indx),' does not exist.']);
            end
            n = obj.M*obj.N/(obj.Kv*obj.Kh);
            if port_indx >= n
                port_indx = mod(port_indx,n);
            end
            pv = obj.M/obj.Kv;
%             ph = obj.N/obj.Kh;
            indx = [mod(port_indx,pv),floor(port_indx/pv)];
            indr = (indx(1,1)*obj.Kv+1):(indx(1,1)+1)*obj.Kv;
            indc = (indx(1,2)*obj.Kh+1):(indx(1,2)+1)*obj.Kh;
        end
        function [indr, indc] = get_element_pos(obj,element_indx)
            if element_indx<0 && element_indx>=obj.num_element
                error(['the port ',num2str(element_indx),' does not exist.']);
            end
            n = obj.M*obj.N;
            if element_indx >= n
                element_indx = mod(element_indx,n);
            end
            pv = obj.M;
%             ph = obj.N/obj.Kh;
            indx = [mod(element_indx,pv),floor(element_indx/pv)];
            indr = (indx(1,1)+1);
            indc = (indx(1,2)+1);
        end
        
        %% get and set function
        function out = get.alpha(obj)
            out = obj.attachedArray.alpha;            
        end
        function out = get.beta(obj)
            out = obj.attachedArray.beta;            
        end
        function out = get.gamma(obj)
            out = obj.attachedArray.gamma;            
        end
        function out = get.attachedType(obj)
            out = obj.attachedArray.attachedType;            
        end
        function out = get.num_element(obj)
            out = obj.M*obj.N*obj.P;
        end
        function out = get.w_m(obj)
            m   = 1:obj.Kv;
            out = 1/sqrt(obj.Kv).*exp(-1j*2*pi*(m'-1).*obj.d_V.*cosd(obj.ele_downtilt));
            if obj.Kv ~= 1
                out = repmat(out,1,obj.Kh);
            end
        end
        function out = get.w_n(obj)
            n   = 1:obj.Kh;
            out = 1/sqrt(obj.Kh).*exp(-1j*2*pi*(n-1).*obj.d_H.*sind(obj.ele_panning));
            if obj.Kh ~= 1
                out = repmat(out,obj.Kv,1);
            end
        end
        function out = get.num_port(obj)
            out = obj.M*obj.N/(obj.Kv*obj.Kh)*obj.P;
        end
        function out = get.pos_element_LCS(obj)
            % get the position of the antenna element in LCS
            [m, n]     = meshgrid(0:obj.M-1,0:obj.N-1);
            out        = zeros(obj.M,obj.N,3);
            out(:,:,2) = (n'-(obj.N-1)/2)*obj.d_H;
            out(:,:,3) = (m'-(obj.M-1)/2)*obj.d_V;
        end
        function out = get.pos_LCS(obj)
            % get the position of the antenna element in LCS
            out = squeeze(obj.attachedArray.pos_panel_LCS(obj.ID(1),obj.ID(2),:));
        end
    end
end