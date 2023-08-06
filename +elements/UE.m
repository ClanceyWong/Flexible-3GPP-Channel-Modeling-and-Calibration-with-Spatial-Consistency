classdef UE < handle
    % UE 
    properties
        ID 
        rand_los
        rand_din
        bIndoor = false
        pos  
        v               % 3 km/h
        theta_v 
        phi_v
        n_fl 
        h_UT 
        antenna 
        link_list
        best_link_id
        couplingloss
        SINR_geometry        % Geometry SINR
        SIR_geometry
        RSRP
        attachedSector = []
        attachBS_pos
        O2IPL
        O2Isigma
        carPL
        fcin
        initPos
    end
    properties(Dependent)
        pos3D
    end
    
    methods
        function obj = UE(varargin)
            if isempty(varargin)
                obj.v       = 3/3.6;   % 3 km/h
                obj.theta_v = 90;
                obj.phi_v   = rand*360-180;
                obj.bIndoor = false;
                obj.n_fl    = 1;
                obj.h_UT    = 1.5;
            else
                % we can provide the parameters;
            end
        end
        
        function [RSRP, couplingloss, pathloss] = RSRP_calc(obj,sector,link,port0)
            tx_power    = sector.attached_BS.Tx_power;
            lambda      = 3e8/sector.frequency;
%             link_params = sector.attached_BS.link_params;
%             link_params = cell2mat(link_params);
%             id          = [link_params(:).UEID];
%             link        = link_params(id == obj.ID);
            KR          = sum(10.^(link.K*0.1));
            [N,M]       = size(link.theta_n_m_ZOD);
            
            theta_n_m_ZOA = link.theta_n_m_ZOA;
            theta_n_m_ZOD = link.theta_n_m_ZOD;
            phi_n_m_AOA   = link.phi_n_m_AOA;
            phi_n_m_AOD   = link.phi_n_m_AOD;
            theta_LOS_ZOA = link.theta_LOS_ZOA;
            theta_LOS_ZOD = link.theta_LOS_ZOD;
            phi_LOS_AOA   = link.phi_LOS_AOA;
            phi_LOS_AOD   = link.phi_LOS_AOD;
            XPR_n_m_      = link.XPR_n_m(:,:,1);
            PHI_n_m       = (2*rand([2, 2, N, M])-1)*pi;
            
            if ~isempty(link.oxygenAbsorp)
                tau_n        = link.tau_n_LOS';
                [~,loss_oxy] = link.oxygenAbsorp.loss(link.d_3D, tau_n);
            else
                loss_oxy = 1;
            end
            
            sector.PHI_n_m= PHI_n_m;
            
            E = obj.antenna.panel(1).num_element/obj.antenna.panel.P;
            U = obj.antenna.num_panel*E;
            
            H     = zeros(U,obj.antenna.panel(1).P);
            H_LOS = zeros(U,obj.antenna.panel(1).P);
            
            % UE params
            r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA).*cosd(phi_n_m_AOA); % the spherical unit vector
            r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA).*sind(phi_n_m_AOA);
            r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA);
            % BS params
            w_m  = sector.antenna.panel(1).w_m.*sector.antenna.panel(1).w_n;
            w_m = w_m(:);
            r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD).*cosd(phi_n_m_AOD); 
            r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD).*sind(phi_n_m_AOD);
            r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD);
            
            XPR_nm(1,1,:,:)  = XPR_n_m_;
            MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m(1,1,:,:));  
            MAT_PHI(1,2,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(1, 2,:,:));
            MAT_PHI(2,1,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(2, 1,:,:));      
            MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m(2, 2,:,:));
            Pwtsqrt          = sqrt(link.loss_blockage).*repmat(sqrt(loss_oxy.*link.Pn.'/M/(KR+1)),1,M);
            
            pos_panelUE = reshape(permute(obj.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            pos_panelBS = reshape(sector.antenna.pos_panel_LCS(1,1,:),3,[]);
            [indr, indc] = sector.antenna.panel(1).get_port_pos(port0);
            for u = 1:U
                panarU = floor((u-1)/E)+1;
                rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,mod((u-1),E)+1))*lambda;
%                 d_rx_s = obj.antenna.R*(pos_panelUE(:,ceil(u/E)) + reshape(permute(obj.antenna.panel(1).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                F_rx_theta_phi_tmp1 = obj.antenna.panel(1).element_list(1).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
                
                d_tx_s             = sector.antenna.R*(pos_panelBS + reshape(permute(sector.antenna.panel(1).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                F_tx_theta_phi_tmp = sector.antenna.panel(1).element_list(1).field_pattern(phi_n_m_AOD, theta_n_m_ZOD);
                F_tx_theta_phi     =  multiprod(multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m,[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
                h_tmp1             = (abs(Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])))).^2;
                H(u,1)             = sum(sum(h_tmp1,2),1);
                
                if obj.antenna.panel.P == 2
                    F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                    F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
                    h_tmp2              = (abs(Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])))).^2;
                    H(u,2)              = sum(sum(h_tmp2,2),1);
                end
            end
            
            if ~link.O2I && link.bLOS
                PHI_LOS  = 2*pi/lambda*sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
                r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
                            sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
                            cosd(theta_LOS_ZOA)];
                r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
                            sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
                            cosd(theta_LOS_ZOD)];
                for u = 1:U
                    panarU = floor((u-1)/E)+1;
                    rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                    d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,mod((u-1),E)+1))*lambda;
%                     d_rx_s = obj.antenna.R*(pos_panelUE(:,ceil(u/E)) + reshape(permute(obj.antenna.panel(1).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                    F_rx_theta_phi_tmp1 = obj.antenna.panel(1).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                    F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                    
                    d_tx_s             = sector.antenna.R*(pos_panelBS + reshape(permute(sector.antenna.panel(1).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                    F_tx_theta_phi_tmp = sector.antenna.panel(1).element_list(1).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
                    F_tx_theta_phi     = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                    H_LOS(u,1)         = (abs(sqrt(KR/(KR+1))*(F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi))).^2; 
                    
                    if obj.antenna.panel.P == 2
                        F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                        F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                        H_LOS(u,2)          = (abs(sqrt(KR/(KR+1))*(F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi))).^2; 
                    
                    end
                end
            end
            RSRP         = -link.PL + link.SF + 10*log10(sum(sum(H_LOS + H))) - 10*log10(U*obj.antenna.panel(1).P) + tx_power;
            couplingloss = -link.PL + link.SF + 10*log10(sum(sum(H_LOS + H))) - 10*log10(U*obj.antenna.panel(1).P);
            pathloss     = -link.PL + link.SF;
        end
        
        function [RSRP, couplingloss] = RSRP_calc_LBW(obj,sector,link,port0)
            tx_power    = sector.attached_BS.Tx_power;
            lambda      = 3e8/sector.frequency;
            if ~isempty(obj.fcin)
                lambdaMT    = 3e8/obj.fcin;
            end
%             link_params = sector.attached_BS.link_params;
%             link_params = cell2mat(link_params);
%             id          = [link_params(:).UEID];
%             link        = link_params(id == obj.ID);
            KR          = sum(10.^(link.K*0.1));
            [N,M]       = size(link.theta_n_m_ZOD);
            
            theta_n_m_ZOA = link.theta_n_m_ZOA;
            theta_n_m_ZOD = link.theta_n_m_ZOD;
            phi_n_m_AOA   = link.phi_n_m_AOA;
            phi_n_m_AOD   = link.phi_n_m_AOD;
            theta_LOS_ZOA = link.theta_LOS_ZOA;
            theta_LOS_ZOD = link.theta_LOS_ZOD;
            phi_LOS_AOA   = link.phi_LOS_AOA;
            phi_LOS_AOD   = link.phi_LOS_AOD;
            XPR_n_m_      = link.XPR_n_m(:,:,1);
            PHI_n_m       = (2*rand([2, 2, N, M])-1)*pi;
            
            if ~isempty(link.oxygenAbsorp)
                tau_n        = link.tau_n_LOS';
                [~,loss_oxy] = link.oxygenAbsorp.loss(link.d_3D, tau_n);
            else
                loss_oxy = 1;
            end
            
            sector.PHI_n_m= PHI_n_m;
            
            E = obj.antenna.panel(1).num_element/obj.antenna.panel.P;
            U = obj.antenna.num_panel*E;
            
            H     = zeros(U,obj.antenna.panel(1).P);
            H_LOS = zeros(U,obj.antenna.panel(1).P);
            
            % UE params
            r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA).*cosd(phi_n_m_AOA); % the spherical unit vector
            r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA).*sind(phi_n_m_AOA);
            r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA);
            % BS params
            w_m  = sector.antenna.panel(1).w_m.*sector.antenna.panel(1).w_n;
            w_m = w_m(:);
            r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD).*cosd(phi_n_m_AOD); 
            r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD).*sind(phi_n_m_AOD);
            r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD);
            
            XPR_nm(1,1,:,:)  = XPR_n_m_;
            MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m(1,1,:,:));  
            MAT_PHI(1,2,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(1, 2,:,:));
            MAT_PHI(2,1,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(2, 1,:,:));      
            MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m(2, 2,:,:));
            Pwtsqrt          = repmat(sqrt(loss_oxy.*link.Pn.'/M/(KR+1)),1,M);
            
            pos_panelUE = reshape(permute(obj.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            pos_panelBS = reshape(sector.antenna.pos_panel_LCS(1,1,:),3,[]);
            [indr, indc] = sector.antenna.panel(1).get_port_pos(port0);
            for u = 1:U
                panarU = floor((u-1)/E)+1;
                rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,mod((u-1),E)+1))*lambda;
%                 d_rx_s = obj.antenna.R*(pos_panelUE(:,ceil(u/E)) + reshape(permute(obj.antenna.panel(1).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                F_rx_theta_phi_tmp1 = obj.antenna.panel(1).element_list(1).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                if exist('lambdaMT','var')
                    F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambdaMT*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                else
                    F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                end
                    F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
                
                d_tx_s             = sector.antenna.R*(pos_panelBS + reshape(permute(sector.antenna.panel(1).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                F_tx_theta_phi_tmp = sector.antenna.panel(1).element_list(1).field_pattern(phi_n_m_AOD, theta_n_m_ZOD);
                if exist('lambdaMT','var')
                    F_tx_theta_phi     =  multiprod(multiprod(exp(1j*2*pi/lambdaMT*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m,[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
                else
                    F_tx_theta_phi     =  multiprod(multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m,[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
                end
                    h_tmp1             = (abs(Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])))).^2;
                
                if obj.antenna.panel.P == 2
                    F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                    if exist('lambdaMT','var')
                        F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambdaMT*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    else
                        F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    end
                    F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
                    h_tmp2              = (abs(Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])))).^2;
                    H(u,2)              = sum(sum(h_tmp2,2),1);
                end
                H(u,1) = sum(sum(h_tmp1,2),1);
            end
            
            if ~link.O2I && link.bLOS
                PHI_LOS  = 2*pi/lambda*sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
                r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
                            sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
                            cosd(theta_LOS_ZOA)];
                r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
                            sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
                            cosd(theta_LOS_ZOD)];
                for u = 1:U
                    panarU = floor((u-1)/E)+1;
                    rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                    d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,mod((u-1),E)+1))*lambda;
%                     d_rx_s = obj.antenna.R*(pos_panelUE(:,ceil(u/E)) + reshape(permute(obj.antenna.panel(1).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                    F_rx_theta_phi_tmp1 = obj.antenna.panel(1).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                    if exist('lambdaMT','var')
                        F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambdaMT*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                    else
                        F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                    end
                    
                    d_tx_s             = sector.antenna.R*(pos_panelBS + reshape(permute(sector.antenna.panel(1).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                    F_tx_theta_phi_tmp = sector.antenna.panel(1).element_list(1).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
                    if exist('lambdaMT','var')
                        F_tx_theta_phi     = (exp(1j*2*pi/lambdaMT*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                    else
                        F_tx_theta_phi     = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                    end
                    H_LOS(u,1)         = (abs(sqrt(KR/(KR+1))*(F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi))).^2; 
                    
                    if obj.antenna.panel.P == 2
                        F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                        if exist('lambdaMT','var')
                            F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambdaMT*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                        else
                            F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                        end
                        H_LOS(u,2)          = (abs(sqrt(KR/(KR+1))*(F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi))).^2; 
                    
                    end
                end
            end
            RSRP         = -link.PL + link.SF + 10*log10(sum(sum(H_LOS + H))) - 10*log10(U*obj.antenna.panel(1).P) + tx_power;
            couplingloss = -link.PL + link.SF + 10*log10(sum(sum(H_LOS + H))) - 10*log10(U*obj.antenna.panel(1).P);
            
        end
        
        function [tau,channel_coef,channel_sum,singular] = get_channel_LBW(obj,sector,t)
            fc          = sector.frequency;
            lambda      = 3e8/fc;
            if ~isempty(obj.fcin)
                lambdaMT    = 3e8/obj.fcin;
            end
            link_params = sector.link_params;
            link_params = cell2mat(link_params);
            id          = [link_params(:).UEID];
            link        = link_params(id == obj.ID);
            KR          = sum(10.^(link.K*0.1));
            [~,M]       = size(link.theta_n_m_ZOD);
            
            theta_n_m_ZOA = link.theta_n_m_ZOA;
            theta_n_m_ZOD = link.theta_n_m_ZOD;
            phi_n_m_AOA   = link.phi_n_m_AOA;
            phi_n_m_AOD   = link.phi_n_m_AOD;
            theta_LOS_ZOA = link.theta_LOS_ZOA;
            theta_LOS_ZOD = link.theta_LOS_ZOD;
            phi_LOS_AOA   = link.phi_LOS_AOA;
            phi_LOS_AOD   = link.phi_LOS_AOD;
            XPR_n_m_      = link.XPR_n_m(:,:,1);
            PHI_n_m       = link.PHI_n_m;
            tau_n         = link.tau_n_LOS;
            indx          = 1:length(tau_n);
            strong_cluster_id = link.strong_cluster_id;
            map_delay         = link.map_delay;         
            ds                = unique(map_delay);
            map = zeros(length(ds),length(map_delay));
            for d = 1:length(ds)
                map(d,:) = map_delay == ds(d);
            end
            N_new             = length(tau_n)+(length(ds)-1)*length(strong_cluster_id);
            tau1              = tau_n(indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)); 
            tau2              = tau_n(strong_cluster_id) + ds';
            tau               = [tau1(:);tau2(:)];
            [tau,si]          = sort(tau);
            
            pu = obj.antenna.panel(1).P;
            E = obj.antenna.panel(1).num_element/pu;
            U = obj.antenna.num_panel*E;
            
            pb = sector.antenna.panel(1).P;
            B = sector.antenna.panel(1).num_port;
            P = sector.antenna.num_panel*B;
            sec = B;
            if sector.antenna.panel(1).P == 2
                sec = B/2;
            end
            
            H     = zeros(U*obj.antenna.panel(1).P,B,N_new);
            H_LOS = zeros(U*obj.antenna.panel(1).P,B,1);
            
            % UE params
            v_vec             = obj.v*[sind(obj.theta_v)*cosd(obj.phi_v);
                                       sind(obj.theta_v)*sind(obj.phi_v);
                                       cosd(obj.theta_v)];
            r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA).*cosd(phi_n_m_AOA); % the spherical unit vector
            r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA).*sind(phi_n_m_AOA);
            r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA);
            phase_vnm         = squeeze(exp(1j*2*pi/lambda*multiprod(r_rx_n_m,v_vec,[1,2],[1,2])*t));
            % BS params
            w_m  = sector.antenna.panel(1).w_m.*sector.antenna.panel(1).w_n;
            w_m = w_m(:);
            r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD).*cosd(phi_n_m_AOD); 
            r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD).*sind(phi_n_m_AOD);
            r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD);
            
            XPR_nm(1,1,:,:)  = XPR_n_m_;
            MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m(1,1,:,:));  
            MAT_PHI(1,2,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(1, 2,:,:));
            MAT_PHI(2,1,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(2, 1,:,:));
            MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m(2, 2,:,:));
            Pwtsqrt          = repmat(sqrt(link.Pn.'/M/(KR+1)),1,M);
            
            pos_panelUE = reshape(permute(obj.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            pos_panelBS = reshape(permute(sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            for u = 0:U-1
                panarU = floor(u/E)+1;
                rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
%                 d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                if exist('lambdaMT','var')
                    F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambdaMT*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                else
                    F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                end
                F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
                if pu == 2
                    F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                    if exist('lambdaMT','var')
                        F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambdaMT*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    else
                        F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    end
                    F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
                end
                
                for p = 0:P-1
                    panar = floor(p/B)+1;
                    port  = mod(p,B);
                    if mod(p,pb)
                        element = 2;
                    else
                        element = 1;
                    end
                    [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
                    d_tx_s             = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                    F_tx_theta_phi_tmp = sector.antenna.panel(panar).element_list(element).field_pattern(phi_n_m_AOD, theta_n_m_ZOD);
                    if exist('lambdaMT','var')
                        F_tx_theta_phi     = multiprod(multiprod(exp(1j*2*pi/lambdaMT*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m,[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
                    else
                        F_tx_theta_phi     = multiprod(multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m,[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
                    end
                    h_tmp1             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])).*phase_vnm;
                    h_tmp11            = sum(h_tmp1((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                    h_tmp12            = h_tmp1(strong_cluster_id,:);
                    h_tmp12            = (h_tmp12*map.').';
                    h_tmp1             = [h_tmp11; h_tmp12(:)];
                    H((u*pu+1),p+1,:)      = h_tmp1(si);
                    if pu == 2
                        h_tmp2             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])).*phase_vnm;
                        h_tmp21            = sum(h_tmp2((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                        h_tmp22            = h_tmp2(strong_cluster_id,:);
                        h_tmp22            = (h_tmp22*map.').';
                        h_tmp2             = [h_tmp21; h_tmp22(:)];
                        H((u+1)*pu,p+1,:)  = h_tmp2(si);
                    end
                end
            end
            
            if ~link.O2I && link.bLOS
                PHI_LOS  = 2*pi/lambda*sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
                r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
                            sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
                            cosd(theta_LOS_ZOA)];
                r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
                            sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
                            cosd(theta_LOS_ZOD)];
                pahse_los = exp(1j*2*pi/lambda*(r_rx_LOS.'*v_vec)*t);
                for u = 0:U-1
                    panarU = floor(u/E)+1;
                    rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                    d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
                    F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                    if exist('lambdaMT','var')
                        F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambdaMT*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                    else
                        F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                    end
                    if pu == 2
                        F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                        if exist('lambdaMT','var')
                            F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambdaMT*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                        else
                            F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                        end
                    end
                    
                    for p = 0:P-1
                        panar = floor(p/B)+1;
                        port  = mod(p,B);
                        if mod(p,pb)
                            element = 2;
                        else
                            element = 1;
                        end
                        [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
                        d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                        F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(element).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
                        if exist('lambdaMT','var')
                            F_tx_theta_phi      = (exp(1j*2*pi/lambdaMT*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                        else
                            F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                        end
                        H_LOS((u*pu+1),p+1,1) = sqrt(KR/(KR+1))*(F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 

                        if pu == 2
                            H_LOS((u+1)*pu,p+1,1) = sqrt(KR/(KR+1))*(F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
                        end
                    end
                    
                end
                H(:,:,1)     = H(:,:,1) + H_LOS;
            end
            channel_coef = H;
            channel_sum  = sum(H,3);
%             singular = 10*log10(svd(channel_sum).^2);
            f = fc; % (sector.BW/10)*((0:9)-(9/2))+
            expf = permute(repmat(exp(-1j*2*pi*tau*f),1,1,size(H,1),size(H,2)),[3,4,1,2]);  % [N, 50]
            channel = repmat(channel_coef,1,1,1,length(f)); % [u,b,N,50]
            channel_freq = permute(sum((channel.*expf),3),[1,2,4,3]);
            channelsum = sum(multiprod(channel_freq,permute(conj(channel_freq),[2,1,3]),[1,2],[1,2]),3);
            singular = 10*log10(svd(channelsum)/numel(f));
        end
        
        function [tau,channel_coef,channel_sum,singular] = get_channel(obj,sector,t,varargin)
            fc          = sector.frequency;
            lambda      = 3e8/fc;
            link_params = sector.link_params;
            link_params = cell2mat(link_params);
            id          = [link_params(:).UEID];
            link        = link_params(id == obj.ID);
            KR          = sum(10.^(link.K*0.1));
            [~,M]       = size(link.theta_n_m_ZOD);
            
            theta_n_m_ZOA = link.theta_n_m_ZOA;
            theta_n_m_ZOD = link.theta_n_m_ZOD;
            phi_n_m_AOA   = link.phi_n_m_AOA;
            phi_n_m_AOD   = link.phi_n_m_AOD;
            theta_LOS_ZOA = link.theta_LOS_ZOA;
            theta_LOS_ZOD = link.theta_LOS_ZOD;
            phi_LOS_AOA   = link.phi_LOS_AOA;
            phi_LOS_AOD   = link.phi_LOS_AOD;
            XPR_n_m_      = link.XPR_n_m(:,:,1);
            PHI_n_m       = link.PHI_n_m;
            tau_n         = link.tau_n_LOS;
            indx          = 1:length(tau_n);
            strong_cluster_id = link.strong_cluster_id;
            map_delay         = link.map_delay;         
            ds                = unique(map_delay);
            map = zeros(length(ds),length(map_delay));
            for d = 1:length(ds)
                map(d,:) = map_delay == ds(d);
            end
            N_new             = length(tau_n)+(length(ds)-1)*length(strong_cluster_id);
            tau1              = tau_n(indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)); 
            tau2              = tau_n(strong_cluster_id) + ds';
            tau               = [tau1(:);tau2(:)];
            [tau,si]          = sort(tau);
            
            pu = obj.antenna.panel(1).P;
            E = obj.antenna.panel(1).num_element/pu;
            U = obj.antenna.num_panel*E;
            
            pb = sector.antenna.panel(1).P;
            B = sector.antenna.panel(1).num_port;
            P = sector.antenna.num_panel*B;
            sec = B;
            if sector.antenna.panel(1).P == 2
                sec = B/2;
            end
            
            H     = zeros(U*obj.antenna.panel(1).P,B,N_new);
            H_LOS = zeros(U*obj.antenna.panel(1).P,B,1);
            
            % UE params
            v_vec             = obj.v*[sind(obj.theta_v)*cosd(obj.phi_v);
                                       sind(obj.theta_v)*sind(obj.phi_v);
                                       cosd(obj.theta_v)];
            r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA).*cosd(phi_n_m_AOA); % the spherical unit vector
            r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA).*sind(phi_n_m_AOA);
            r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA);
            phase_vnm         = squeeze(exp(1j*2*pi/lambda*multiprod(r_rx_n_m,v_vec,[1,2],[1,2])*t));
            % BS params
            w_m  = sector.antenna.panel(1).w_m.*sector.antenna.panel(1).w_n;
            w_m = w_m(:);
            r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD).*cosd(phi_n_m_AOD); 
            r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD).*sind(phi_n_m_AOD);
            r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD);
            
            XPR_nm(1,1,:,:)  = XPR_n_m_;
            MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m(1,1,:,:));  
            MAT_PHI(1,2,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(1, 2,:,:));
            MAT_PHI(2,1,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(2, 1,:,:));
            MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m(2, 2,:,:));
            Pwtsqrt          = repmat(sqrt(link.Pn.'/M/(KR+1)),1,M);
            
            pos_panelUE = reshape(permute(obj.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            pos_panelBS = reshape(permute(sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            for u = 0:U-1
                panarU = floor(u/E)+1;
                rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
%                 d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
                if pu == 2
                    F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                    F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
                end
                
                for p = 0:P-1
                    panar = floor(p/B)+1;
                    port  = mod(p,B);
                    if mod(p,pb)
                        element = 2;
                    else
                        element = 1;
                    end
                    [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
                    d_tx_s             = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                    F_tx_theta_phi_tmp = sector.antenna.panel(panar).element_list(element).field_pattern(phi_n_m_AOD, theta_n_m_ZOD);
                    F_tx_theta_phi     = multiprod(multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m,[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
                    h_tmp1             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])).*phase_vnm;
                    h_tmp11            = sum(h_tmp1((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                    h_tmp12            = h_tmp1(strong_cluster_id,:);
                    h_tmp12            = (h_tmp12*map.').';
                    h_tmp1             = [h_tmp11; h_tmp12(:)];
                    H((u*pu+1),p+1,:)      = h_tmp1(si);
                    if pu == 2
                        h_tmp2             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])).*phase_vnm;
                        h_tmp21            = sum(h_tmp2((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                        h_tmp22            = h_tmp2(strong_cluster_id,:);
                        h_tmp22            = (h_tmp22*map.').';
                        h_tmp2             = [h_tmp21; h_tmp22(:)];
                        H((u+1)*pu,p+1,:)  = h_tmp2(si);
                    end
                end
            end
            
            if ~link.O2I && link.bLOS
                d_LOS = sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
                if ~isempty(varargin) && varargin{1}
                    phi_GR_AOD = phi_LOS_AOD;
                    phi_GR_AOA = phi_GR_AOD + 180;
                    theta_GR_ZOD = 180 - atand((sqrt(sum((obj.pos3D(1:2)- [link.BS_pos_wrap]').^2)))/(obj.pos3D(3)+sector.attached_BS.h_BS));
                    d_GR  = sqrt(sum((obj.pos3D(1:2)- [link.BS_pos_wrap]').^2) + (obj.pos3D(3)+sector.attached_BS.h_BS).^2 );
                    theta_GR_ZOA = theta_GR_ZOD;
                    mGR = addition_components.GroundReflection(9);
                    PHI_GR = mGR.GetPhi(fc,theta_GR_ZOD);
                    r_rx_GR = [sind(theta_GR_ZOA)*cosd(phi_GR_AOA);  % the spherical unit vector
                               sind(theta_GR_ZOA)*sind(phi_GR_AOA); 
                               cosd(theta_GR_ZOA)];
                    r_tx_GR = [sind(theta_GR_ZOD)*cosd(phi_GR_AOD);  % the spherical unit vector
                               sind(theta_GR_ZOD)*sind(phi_GR_AOD); 
                               cosd(theta_GR_ZOD)];
                    pahse_GR = exp(-1j*2*pi*d_GR/lambda)*exp(1j*2*pi/lambda*(r_rx_GR.'*v_vec)*t);
                    
                    for u = 0:U-1
                        panarU = floor(u/E)+1;
                        rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                        d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
                        F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_GR_AOA, theta_GR_ZOA);
                        F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_GR'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                        if pu == 2
                            F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_GR_AOA, theta_GR_ZOA);
                            F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_GR'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                        end

                        for p = 0:P-1
                            panar = floor(p/B)+1;
                            port  = mod(p,B);
                            if mod(p,pb)
                                element = 2;
                            else
                                element = 1;
                            end
                            [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
                            d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                            F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(element).field_pattern(phi_GR_AOD, theta_GR_ZOD);
                            F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_GR'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                            H_LOS((u*pu+1),p+1,1) = (d_LOS/d_GR)*(F_rx_theta_phi1.'*PHI_GR*F_tx_theta_phi).*pahse_GR; 

                            if pu == 2
                                H_LOS((u+1)*pu,p+1,1) = (d_LOS/d_GR)*(F_rx_theta_phi2.'*PHI_GR*F_tx_theta_phi).*pahse_GR; 
                            end
                        end

                    end
                end
                
                PHI_LOS  = 2*pi/lambda*d_LOS;
                r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
                            sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
                            cosd(theta_LOS_ZOA)];
                r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
                            sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
                            cosd(theta_LOS_ZOD)];
                pahse_los = exp(1j*2*pi/lambda*(r_rx_LOS.'*v_vec)*t);
                for u = 0:U-1
                    panarU = floor(u/E)+1;
                    rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                    d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
                    F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                    F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
                    if pu == 2
                        F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                        F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
                    end
                    
                    for p = 0:P-1
                        panar = floor(p/B)+1;
                        port  = mod(p,B);
                        if mod(p,pb)
                            element = 2;
                        else
                            element = 1;
                        end
                        [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
                        d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                        F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(element).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
                        F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
                        H_LOS((u*pu+1),p+1,1) = H_LOS((u*pu+1),p+1,1) + (F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 

                        if pu == 2
                            H_LOS((u+1)*pu,p+1,1) = H_LOS((u+1)*pu,p+1,1) + (F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
                        end
                    end
                    
                end
                H(:,:,1)     = H(:,:,1) + sqrt(KR/(KR+1))*H_LOS;
            end
%             if ~link.O2I && link.bLOS
%                 PHI_LOS  = 2*pi/lambda*sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
%                 r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
%                             sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
%                             cosd(theta_LOS_ZOA)];
%                 r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
%                             sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
%                             cosd(theta_LOS_ZOD)];
%                 pahse_los = exp(1j*2*pi/lambda*(r_rx_LOS.'*v_vec)*t);
%                 for u = 0:U-1
%                     panarU = floor(u/E)+1;
%                     rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
%                     d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
%                     F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
%                     F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
%                     if pu == 2
%                         F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
%                         F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
%                     end
%                     
%                     for p = 0:P-1
%                         panar = floor(p/B)+1;
%                         port  = mod(p,B);
%                         if mod(p,pb)
%                             element = 2;
%                         else
%                             element = 1;
%                         end
%                         [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
%                         d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
%                         F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(element).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
%                         F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
%                         H_LOS((u*pu+1),p+1,1) = sqrt(KR/(KR+1))*(F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
% 
%                         if pu == 2
%                             H_LOS((u+1)*pu,p+1,1) = sqrt(KR/(KR+1))*(F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
%                         end
%                     end
%                     
%                 end
%                 H(:,:,1)     = H(:,:,1) + H_LOS;
%             end
            channel_coef = H;
            channel_sum  = sum(H,3);
%             singular = 10*log10(svd(channel_sum).^2);
            f = fc; % (sector.BW/10)*((0:9)-(9/2))+
            expf = permute(repmat(exp(-1j*2*pi*(tau)*f),1,1,size(H,1),size(H,2)),[3,4,1,2]);  % [N, 50]
            channel = repmat(channel_coef,1,1,1,length(f)); % [u,b,N,50]
            channel_freq = permute(sum((channel.*expf),3),[1,2,4,3]);
            channelsum = sum(multiprod(channel_freq,permute(conj(channel_freq),[2,1,3]),[1,2],[1,2]),3);
            singular = 10*log10(svd(channelsum)/numel(f));
        end
        
        
        function [tau,channel_coef,channel_sum,singular] = get_channel_v2(obj,sector,t,varargin)
            fc          = sector.frequency;
            lambda      = 3e8/fc;
            link_params = sector.link_params;
            link_params = cell2mat(link_params);
            id          = [link_params(:).UEID];
            link        = link_params(id == obj.ID);
            KR          = sum(10.^(link.K*0.1));
            [~,M]       = size(link.theta_n_m_ZOD);
            
            theta_n_m_ZOA = link.theta_n_m_ZOA;
            theta_n_m_ZOD = link.theta_n_m_ZOD;
            phi_n_m_AOA   = link.phi_n_m_AOA;
            phi_n_m_AOD   = link.phi_n_m_AOD;
            theta_LOS_ZOA = link.theta_LOS_ZOA;
            theta_LOS_ZOD = link.theta_LOS_ZOD;
            phi_LOS_AOA   = link.phi_LOS_AOA;
            phi_LOS_AOD   = link.phi_LOS_AOD;
            XPR_n_m_      = link.XPR_n_m(:,:,1);
            PHI_n_m       = link.PHI_n_m;
            tau_n         = link.tau_n_LOS;
            indx          = 1:length(tau_n);
            strong_cluster_id = link.strong_cluster_id;
            map_delay         = link.map_delay;         
            ds                = unique(map_delay);
            map = zeros(length(ds),length(map_delay));
            for d = 1:length(ds)
                map(d,:) = map_delay == ds(d);
            end
            N_new             = length(tau_n)+(length(ds)-1)*length(strong_cluster_id);
            tau1              = tau_n(indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)); 
            tau2              = tau_n(strong_cluster_id) + ds';
            tau               = [tau1(:);tau2(:)];
            [tau,si]          = sort(tau);
            
            pu = obj.antenna.panel(1).P;
            E = obj.antenna.panel(1).num_element/pu;
            U = obj.antenna.num_panel*E;
            
            pb = sector.antenna.panel(1).P;
            B = sector.antenna.panel(1).num_element/pb;
            P = sector.antenna.num_panel*B;
            sec = B;
            if sector.antenna.panel(1).P == 2
                sec = B/2;
            end
            
            H     = zeros(U*obj.antenna.panel(1).P,P*sector.antenna.panel(1).P,N_new);
            H_LOS = zeros(U*obj.antenna.panel(1).P,P*sector.antenna.panel(1).P,1);
            
            % UE params
            v_vec             = obj.v*[sind(obj.theta_v)*cosd(obj.phi_v);
                                       sind(obj.theta_v)*sind(obj.phi_v);
                                       cosd(obj.theta_v)];
            r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA).*cosd(phi_n_m_AOA); % the spherical unit vector
            r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA).*sind(phi_n_m_AOA);
            r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA);
            phase_vnm         = squeeze(exp(1j*2*pi/lambda*multiprod(r_rx_n_m,v_vec,[1,2],[1,2])*t));
            % BS params
            w_m  = sector.antenna.panel(1).w_m.*sector.antenna.panel(1).w_n;
            w_m = w_m(:);
            r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD).*cosd(phi_n_m_AOD); 
            r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD).*sind(phi_n_m_AOD);
            r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD);
            
            XPR_nm(1,1,:,:)  = XPR_n_m_;
            MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m(1,1,:,:));  
            MAT_PHI(1,2,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(1, 2,:,:));
            MAT_PHI(2,1,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m(2, 1,:,:));
            MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m(2, 2,:,:));
            Pwtsqrt          = repmat(sqrt(link.Pn.'/M/(KR+1)),1,M);
            
            pos_panelUE = reshape(permute(obj.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            pos_panelBS = reshape(permute(sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);
            for u = 0:U-1
                panarU = floor(u/E)+1;
                rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
%                 d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]))*lambda;
                F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                F_rx_theta_phi_tmp1 =  multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp1,[1 2],[1 2]);
                F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
                if pu == 2
                    F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_n_m_AOA, theta_n_m_ZOA);
                    F_rx_theta_phi_tmp2 = multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),F_rx_theta_phi_tmp2,[1 2],[1 2]);
                    F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
                end
                
                for p = 0:P-1
                    panar = floor(p/B)+1;
                    port  = mod(p,B);
                    [indr, indc] = sector.antenna.panel(panar).get_element_pos(port);
                    d_tx_s             = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                    F_tx_theta_phi_tmp = sector.antenna.panel(panar).element_list(1).field_pattern(phi_n_m_AOD, theta_n_m_ZOD);
                    F_tx_theta_phi     = multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),F_tx_theta_phi_tmp,[1 2],[1 2]);
                    h_tmp1             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])).*phase_vnm;
                    h_tmp11            = sum(h_tmp1((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                    h_tmp12            = h_tmp1(strong_cluster_id,:);
                    h_tmp12            = (h_tmp12*map.').';
                    h_tmp1             = [h_tmp11; h_tmp12(:)];
                    H((u*pu+1),p*pb+1,:)      = h_tmp1(si);
                    if pu == 2
                        h_tmp2             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])).*phase_vnm;
                        h_tmp21            = sum(h_tmp2((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                        h_tmp22            = h_tmp2(strong_cluster_id,:);
                        h_tmp22            = (h_tmp22*map.').';
                        h_tmp2             = [h_tmp21; h_tmp22(:)];
                        H((u+1)*pu,p*pb+1,:)  = h_tmp2(si);
                    end
                    
                    if pb == 2
                        F_tx_theta_phi_tmp = sector.antenna.panel(panar).element_list(2).field_pattern(phi_n_m_AOD, theta_n_m_ZOD);
                        F_tx_theta_phi2     = multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),F_tx_theta_phi_tmp,[1 2],[1 2]);
                        h_tmp3             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi2,[1 2],[1 2])).*phase_vnm;
                        h_tmp31            = sum(h_tmp3((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                        h_tmp32            = h_tmp3(strong_cluster_id,:);
                        h_tmp32            = (h_tmp32*map.').';
                        h_tmp3             = [h_tmp31; h_tmp32(:)];
                        H((u*pu+1),(p+1)*pb,:)      = h_tmp3(si);
                        if pu == 2
                            h_tmp4             = Pwtsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi2,[1 2],[1 2])).*phase_vnm;
                            h_tmp41            = sum(h_tmp4((indx~=strong_cluster_id(1)&indx~=strong_cluster_id(2)),:),2);
                            h_tmp42            = h_tmp4(strong_cluster_id,:);
                            h_tmp42            = (h_tmp42*map.').';
                            h_tmp4             = [h_tmp41; h_tmp42(:)];
                            H((u+1)*pu,(p+1)*pb,:)  = h_tmp4(si);
                        end
                    end
                end
            end
            
            if ~link.O2I && link.bLOS
                d_LOS = sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
                if ~isempty(varargin) && varargin{1}
                    phi_GR_AOD = phi_LOS_AOD;
                    phi_GR_AOA = phi_GR_AOD + 180;
                    theta_GR_ZOD = 180 - atand((sqrt(sum((obj.pos3D(1:2)- [link.BS_pos_wrap]').^2)))/(obj.pos3D(3)+sector.attached_BS.h_BS));
                    d_GR  = sqrt(sum((obj.pos3D(1:2)- [link.BS_pos_wrap]').^2) + (obj.pos3D(3)+sector.attached_BS.h_BS).^2 );
                    theta_GR_ZOA = theta_GR_ZOD;
                    mGR = addition_components.GroundReflection(9);
                    PHI_GR = mGR.GetPhi(fc,theta_GR_ZOD);
                    r_rx_GR = [sind(theta_GR_ZOA)*cosd(phi_GR_AOA);  % the spherical unit vector
                               sind(theta_GR_ZOA)*sind(phi_GR_AOA); 
                               cosd(theta_GR_ZOA)];
                    r_tx_GR = [sind(theta_GR_ZOD)*cosd(phi_GR_AOD);  % the spherical unit vector
                               sind(theta_GR_ZOD)*sind(phi_GR_AOD); 
                               cosd(theta_GR_ZOD)];
                    pahse_GR = exp(-1j*2*pi*d_GR/lambda)*exp(1j*2*pi/lambda*(r_rx_GR.'*v_vec)*t);
                    
                    for u = 0:U-1
                        panarU = floor(u/E)+1;
                        rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                        d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
                        F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_GR_AOA, theta_GR_ZOA);
                        F_rx_theta_phi1     =  (exp(1j*2*pi/lambda*(r_rx_GR'*d_rx_s)))*F_rx_theta_phi_tmp1;
                        if pu == 2
                            F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_GR_AOA, theta_GR_ZOA);
                            F_rx_theta_phi2     =  (exp(1j*2*pi/lambda*(r_rx_GR'*d_rx_s)))*F_rx_theta_phi_tmp2;
                        end

                        for p = 0:P-1
                            panar = floor(p/B)+1;
                            port  = mod(p,B);
                            [indr, indc] = sector.antenna.panel(panar).get_element_pos(port);
                            d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                            F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(1).field_pattern(phi_GR_AOD, theta_GR_ZOD);
                            F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_GR'*d_tx_s)))*F_tx_theta_phi_tmp;
                            H_LOS((u*pu+1),p*pb+1,1) = (d_LOS/d_GR)*(F_rx_theta_phi1.'*PHI_GR*F_tx_theta_phi).*pahse_GR; 

                            if pu == 2
                                H_LOS((u+1)*pu,p*pb+1,1) = (d_LOS/d_GR)*(F_rx_theta_phi2.'*PHI_GR*F_tx_theta_phi).*pahse_GR; 
                            end
                            if pb == 2
                                F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(2).field_pattern(phi_GR_AOD, theta_GR_ZOD);
                                F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_GR'*d_tx_s)))*F_tx_theta_phi_tmp;
                                H_LOS((u*pu+1),(p+1)*pb,1) = (d_LOS/d_GR)*(F_rx_theta_phi1.'*PHI_GR*F_tx_theta_phi).*pahse_GR; 

                                if pu == 2
                                    H_LOS((u+1)*pu,(p+1)*pb,1) = (d_LOS/d_GR)*(F_rx_theta_phi2.'*PHI_GR*F_tx_theta_phi).*pahse_GR; 
                                end
                            end
                        end

                    end
                end
                
                PHI_LOS  = 2*pi/lambda*d_LOS;
                r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
                            sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
                            cosd(theta_LOS_ZOA)];
                r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
                            sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
                            cosd(theta_LOS_ZOD)];
                pahse_los = exp(1j*2*pi/lambda*(r_rx_LOS.'*v_vec)*t);
                for u = 0:U-1
                    panarU = floor(u/E)+1;
                    rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
                    d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
                    F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                    F_rx_theta_phi1     =  (exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)))*F_rx_theta_phi_tmp1;
                    if pu == 2
                        F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
                        F_rx_theta_phi2     =  (exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)))*F_rx_theta_phi_tmp2;
                    end
                    
                    for p = 0:P-1
                        panar = floor(p/B)+1;
                        port  = mod(p,B);
                        [indr, indc] = sector.antenna.panel(panar).get_element_pos(port);
                        d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
                        F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(1).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
                        F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s)))*F_tx_theta_phi_tmp;
                        H_LOS((u*pu+1),p*pb+1,1) = H_LOS((u*pu+1),p*pb+1,1) + (F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 

                        if pu == 2
                            H_LOS((u+1)*pu,p*pb+1,1) = H_LOS((u+1)*pu,p*pb+1,1) + (F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
                        end
                        
                        if pb == 2
                            F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(2).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
                            F_tx_theta_phi2      = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s)))*F_tx_theta_phi_tmp;
                            H_LOS((u*pu+1),(p+1)*pb,1) = H_LOS((u*pu+1),(p+1)*pb,1) + (F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi2).*pahse_los; 

                            if pu == 2
                                H_LOS((u+1)*pu,(p+1)*pb,1) = H_LOS((u+1)*pu,(p+1)*pb,1) + (F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi2).*pahse_los; 
                            end
                        end
                    end
                    
                end
                H(:,:,1)     = H(:,:,1) + sqrt(KR/(KR+1))*H_LOS;
            end
%             if ~link.O2I && link.bLOS
%                 PHI_LOS  = 2*pi/lambda*sqrt(sum((obj.pos3D- [link.BS_pos_wrap, sector.attached_BS.h_BS]').^2));
%                 r_rx_LOS = [sind(theta_LOS_ZOA)*cosd(phi_LOS_AOA);  % the spherical unit vector
%                             sind(theta_LOS_ZOA)*sind(phi_LOS_AOA); 
%                             cosd(theta_LOS_ZOA)];
%                 r_tx_LOS = [sind(theta_LOS_ZOD)*cosd(phi_LOS_AOD);  % the spherical unit vector
%                             sind(theta_LOS_ZOD)*sind(phi_LOS_AOD); 
%                             cosd(theta_LOS_ZOD)];
%                 pahse_los = exp(1j*2*pi/lambda*(r_rx_LOS.'*v_vec)*t);
%                 for u = 0:U-1
%                     panarU = floor(u/E)+1;
%                     rx_pos = reshape(permute(obj.antenna.panel(panarU).pos_element_LCS,[3,1,2]),3,[]);
%                     d_rx_s = obj.antenna.R*(pos_panelUE(:,panarU) + rx_pos(:,u+1))*lambda;
%                     F_rx_theta_phi_tmp1 = obj.antenna.panel(panarU).element_list(1).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
%                     F_rx_theta_phi1     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp1;
%                     if pu == 2
%                         F_rx_theta_phi_tmp2 = obj.antenna.panel(panarU).element_list(2).field_pattern(phi_LOS_AOA, theta_LOS_ZOA);
%                         F_rx_theta_phi2     =  sum(exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s)),2)*F_rx_theta_phi_tmp2;
%                     end
%                     
%                     for p = 0:P-1
%                         panar = floor(p/B)+1;
%                         port  = mod(p,B);
%                         if mod(p,pb)
%                             element = 2;
%                         else
%                             element = 1;
%                         end
%                         [indr, indc] = sector.antenna.panel(panar).get_port_pos(port);
%                         d_tx_s              = sector.antenna.R*(pos_panelBS(:,panar) + reshape(permute(sector.antenna.panel(panar).pos_element_LCS(indr,indc,:),[3,1,2]),3,[]))*lambda;
%                         F_tx_theta_phi_tmp  = sector.antenna.panel(panar).element_list(element).field_pattern(phi_LOS_AOD, theta_LOS_ZOD);
%                         F_tx_theta_phi      = (exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))*w_m)*F_tx_theta_phi_tmp;
%                         H_LOS((u*pu+1),p+1,1) = sqrt(KR/(KR+1))*(F_rx_theta_phi1.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
% 
%                         if pu == 2
%                             H_LOS((u+1)*pu,p+1,1) = sqrt(KR/(KR+1))*(F_rx_theta_phi2.'*[exp(-1j*PHI_LOS), 0; 0, -exp(-1j*PHI_LOS)]*F_tx_theta_phi).*pahse_los; 
%                         end
%                     end
%                     
%                 end
%                 H(:,:,1)     = H(:,:,1) + H_LOS;
%             end
            channel_coef = H;
            channel_sum  = sum(H,3);
%             singular = 10*log10(svd(channel_sum).^2);
            f = fc; % (sector.BW/10)*((0:9)-(9/2))+
            expf = permute(repmat(exp(-1j*2*pi*(tau)*f),1,1,size(H,1),size(H,2)),[3,4,1,2]);  % [N, 50]
            channel = repmat(channel_coef,1,1,1,length(f)); % [u,b,N,50]
            channel_freq = permute(sum((channel.*expf),3),[1,2,4,3]);
            channelsum = sum(multiprod(channel_freq,permute(conj(channel_freq),[2,1,3]),[1,2],[1,2]),3);
            singular = 10*log10(svd(channelsum)/numel(f));
        end
        
        function value = get.pos3D(obj)
            value = [obj.pos';obj.h_UT];
        end
        
        function plot_link(obj)
            if isempty(obj.attachedSector)
                % do nothing
            else
                BSpos = [obj.attachBS_pos, obj.attachedSector.attached_BS.h_BS].';
                UEpos = obj.pos3D;
                figure(1); hold on;
                plot3([UEpos(1),BSpos(1)], [UEpos(2),BSpos(2)], [UEpos(3),BSpos(3)],'color',obj.attachedSector.lineColor);hold on;
            end
       end
    end
end