classdef Link38900_e00 < handle
    % LINK
    properties
        fastfading_enable = false
        sector
        BS_pos_wrap
        UE
        O2I = false
        scenario
        fc
        d_2D
        d_2D_in
        d_2D_out
        d_3D
        
        phi_LOS_AOD
        theta_LOS_ZOD
        phi_LOS_AOA
        theta_LOS_ZOA
        bLOS
        
        PL
        sigma_SF
        CouplingLoss
        
        SF
        K
        DS 
        ASD 
        ASA 
        mu_lgZSD
        sigma_lgZSD
        mu_offset_ZOD
        ZSD 
        ZSA
        r_tau 
        mu_XPR 
        sigma_XPR 
        N        % Number of clusters 
        M        % Number of rays per cluster
        zeta     % Per cluster shadowing std
        c_DS
        c_ASA 
        c_ZSA 
        c_ASD
        
        tau_order
        tau_n
        tau_n_LOS
        Pn 
        Pn_LOS
        keep
        N_new
        strong_cluster_id
        map_delay
        
        phi_n_m_AOA
        phi_n_m_AOD
        theta_n_m_ZOA
        theta_n_m_ZOD
        XPR_n_m
        PHI_n_m
        PHI_LOS 
        
        oxygenAbsorp
        
        spatialConsist_indx
        
        loss_blockage
    end
    
    methods
        function obj = Link38900_e00(BSsector, UE, scenario, fastfading_enable, tab, varargin)
            obj.sector = BSsector;
            obj.UE = UE;
            obj.scenario = scenario;
            obj.fc = BSsector.frequency;
            obj.fastfading_enable = fastfading_enable;
            if UE.bIndoor
                obj.O2I = true;
            end
            
            if strcmp(scenario.name,'3D-RMa') ||strcmp(scenario.name,'3D-UMa')||strcmp(scenario.name,'3D-UMi')
                [~, indx] = min(sum(((repmat(UE.pos, size(BSsector.Pos_wrapped,1), 1)-BSsector.Pos_wrapped).^2),2));
                obj.BS_pos_wrap = BSsector.Pos_wrapped(indx,:);
            elseif strcmp(scenario.name,'3D-InH')
                obj.BS_pos_wrap = BSsector.Position;
            end
            
            
           %% step 1
            x = UE.pos(1)-obj.BS_pos_wrap(1);
            y = UE.pos(2)-obj.BS_pos_wrap(2);
            z = UE.h_UT-BSsector.h_BS;
            obj.d_2D = sqrt(x.^2+y.^2);
            obj.d_3D = sqrt(x.^2+y.^2+z.^2);
%             obj.d_2D_in = min(UE.rand_din(1),obj.d_2D);
%             obj.d_2D_in = min(UE.rand_din(BSsector.attached_BS.ID),obj.d_2D);
%             if obj.O2I
            if strcmp(scenario.name,'3D-RMa') ||strcmp(scenario.name,'3D-UMa')||strcmp(scenario.name,'3D-UMi')
                if length(UE.rand_din) == 1
                    obj.d_2D_in  = UE.rand_din(1);
                else
                    obj.d_2D_in  = UE.rand_din(BSsector.ID);
                end
                obj.d_2D_out = obj.d_2D - obj.d_2D_in;
            elseif strcmp(scenario.name,'3D-InH')
                obj.d_2D_in  = obj.d_2D;
                obj.d_2D_out = 0;
            end
%             else
%                 obj.d_2D_out = obj.d_2D;
%             end
            [obj.phi_LOS_AOD, obj.theta_LOS_ZOD] = cart2sph_(x, y, z);
            [obj.phi_LOS_AOA, obj.theta_LOS_ZOA] = cart2sph_(-x, -y, -z);
            
           %% step 2
            PrLOS = obj.los_probability();
            if ~isempty(obj.sector.SpatialConsistency)
                randlos = obj.sector.SpatialConsistency.LOSstate.rand([obj.UE.pos,obj.UE.h_UT]');
            else
                randlos = UE.rand_los(BSsector.ID);
            end
            if randlos < PrLOS
                obj.bLOS = true;
                obj.spatialConsist_indx = 1;
            else
                obj.bLOS = false;
                obj.spatialConsist_indx = 2;
            end
            if obj.O2I
                obj.spatialConsist_indx = 3;
            end
            
           %% step 3
            % obj.bLOS = true;
            obj.pathloss();
            
           %% step 4
            obj.large_scale_para(tab);
            obj.ASA = min(obj.ASA, 104);
            obj.ASD = min(obj.ASD, 104);
            obj.ZSA = min(obj.ZSA, 52);
            obj.ZSD = min(obj.ZSD, 52);
            
            % the following steps is needed only when fast fading is enable
            if fastfading_enable
              %% step 5: Generate cluster delays
                if ~isempty(obj.sector.SpatialConsistency)
                    Xn = obj.sector.SpatialConsistency.delay{obj.spatialConsist_indx}.rand([obj.UE.pos,obj.UE.h_UT]');
                else
                    Xn = rand(obj.N,1);
                end
                tau_n_tmp = -obj.r_tau*obj.DS*log(Xn.');
                [obj.tau_n, obj.tau_order] = sort(tau_n_tmp - min(tau_n_tmp));
                obj.tau_n_LOS = obj.tau_n;
                if obj.bLOS && ~obj.O2I
                    C_tau = 0.7705-0.0433*obj.K+0.0002*obj.K^2+0.000017*obj.K^3;
                    obj.tau_n_LOS = obj.tau_n/C_tau;  % not to be used in cluster power generation
                end
                
                if ~isempty(varargin) 
                    if strcmp(varargin{1},'oxygen absorption')
                        if ~obj.O2I && obj.bLOS
                            tau_delta = 0;
                        else
                            tau_delta = min(tau_n_tmp);
                        end
                        obj.oxygenAbsorp = addition_components.Oxygen_Absorption(obj.fc/1e9, tau_delta);
                    end
                end
              %% step 6: Generate cluster powers
                if ~isempty(obj.sector.SpatialConsistency)
                    Zn = obj.sector.SpatialConsistency.shadowing{obj.spatialConsist_indx}.randn([obj.UE.pos,obj.UE.h_UT]')*obj.zeta;
                    Zn = Zn(obj.tau_order);
                else
                    Zn = randn(obj.N,1)*obj.zeta;
                end
                Pn_tmp = exp(-obj.tau_n.*(obj.r_tau-1)./obj.r_tau./obj.DS).*(10.^(-Zn.'/10));
                obj.Pn = Pn_tmp/sum(Pn_tmp);
                obj.Pn_LOS = obj.Pn;
                if obj.bLOS && ~obj.O2I
                    KR = 10^(obj.K/10);
                    P1_LOS = KR/(KR+1);
                    obj.Pn_LOS = (1/(KR+1)).*(Pn_tmp./sum(Pn_tmp));
                    obj.Pn_LOS(1) = obj.Pn_LOS(1) + P1_LOS;  % not be used in equation step 11(7.3-22)
                end
%                 obj.cluster_power_per_ray = obj.Pn/obj.M;
                if ~isempty(obj.sector.SpatialConsistency)
                    keep       = true(1,numel(obj.Pn_LOS));
                else
                    keep       = ((10*log10(obj.Pn_LOS/max(obj.Pn_LOS))) >= -25);
                end
                if sum(keep) == 1
                    keep(obj.Pn_LOS == max(obj.Pn_LOS(~keep))) = true;
                end
                obj.keep  = keep;
                obj.tau_n  = obj.tau_n(keep);
                obj.tau_n_LOS  = obj.tau_n_LOS(keep);
                obj.Pn     = obj.Pn(keep);
                obj.Pn_LOS = obj.Pn_LOS(keep);
                obj.N_new  = numel(obj.Pn);
                [~,tmp]    = sort((1./obj.Pn));
                obj.strong_cluster_id = tmp(1:min(2,numel(tmp)));  % cluster number is 2
                obj.map_delay = zeros(1,obj.M);
                if isempty(obj.sector.SpatialConsistency)
                    if ~isempty(obj.c_DS)
                        obj.map_delay([9,10,11,12,17,18]) = obj.c_DS*1.28e-9;
                        obj.map_delay(13:16) = obj.c_DS*2.56e-9;
                    else
                        obj.map_delay(9:18) = 3.91e-9;
                    end
                end
                    
               %% step 7: Generate arrival angles and departure angles for both azimuth and elevation
                obj.AOA_calc();
                obj.AOD_calc();
                obj.ZOA_calc();
                obj.ZOD_calc();
                
               %% step 8: Coupling of rays within a cluster for both azimuth and elevation
                obj.RandomCouplingRays();
                
              %% step 9: Generate XPRs
                if ~isempty(obj.sector.SpatialConsistency)
                    rn = obj.sector.SpatialConsistency.XPR{obj.spatialConsist_indx}.randn([obj.UE.pos,obj.UE.h_UT]');
                    rn = reshape(rn,obj.N,[]);
                    rn = rn(obj.tau_order,:);
                    rn = rn(obj.keep,:);
                else
                    rn = randn(obj.N_new, obj.M);
                end
                X_n_m = obj.sigma_XPR*rn+obj.mu_XPR;
%                 X_n_m(X_n_m<0) = 0;
                obj.XPR_n_m = 10.^(X_n_m/10);
                
                if ~isempty(obj.sector.blockage)
                    inv_R  = obj.UE.antenna.R;
                    xyz    = sph2cart_(obj.phi_n_m_AOA, obj.theta_n_m_ZOA);
                    tmp    = multiprod(inv_R,xyz,[1 2],1);
                    tmp1   = tmp(1,:,:); 
                    tmp2   = tmp(2,:,:); 
                    tmp3   = tmp(3,:,:);
                    [phiLCS, thetaLCS] = cart2sph_(tmp1, tmp2, tmp3);
                    for ik = 1:size(phiLCS,2)
                        [lossSelf(ik,1), loss(ik,1)] = obj.sector.blockage.attenuation(thetaLCS(1,ik,1),phiLCS(1,ik,1),obj.theta_n_m_ZOA(ik,1),obj.phi_n_m_AOA(ik,1),[obj.UE.pos,obj.UE.h_UT]',obj.spatialConsist_indx);
                    end
                    obj.loss_blockage = repmat(10.^((lossSelf+loss)/10),1,size(obj.phi_n_m_AOA,2));
                else
                    obj.loss_blockage = ones(size(obj.phi_n_m_AOA));
                end
                
                % The outcome of Steps 1-9 shall be identical for all the links from co-sited sectors to a UT.
                
              %% step 10: Draw initial random phases
                if ~isempty(obj.sector.SpatialConsistency)
                    rnp = obj.sector.SpatialConsistency.randomPhase{obj.spatialConsist_indx}.rand([obj.UE.pos,obj.UE.h_UT]');
                    rnp = reshape(rnp,2,2,obj.N,[]);
                    rnp = rnp(:,:,obj.tau_order,:);
                    rnp = rnp(:,:,obj.keep,:);
                    obj.PHI_n_m = (2*rnp-1)*pi;
                else
                    obj.PHI_n_m = (2*rand([2, 2, obj.N_new, obj.M])-1)*pi;
                end
            end
        end
        
        
        % caculator the LOS probability
        % The LOS probability is derived with assuming antenna heights of 
        % 3m for indoor, 10m for UMi, and 25m for UMa
        function Pr_LOS = los_probability(obj)
            d_2Dout = obj.d_2D_out;
            h_UT = obj.UE.h_UT;
            if strcmp(obj.scenario.name,'3D-UMa')
                if d_2Dout <= 18
                    Pr_LOS = 1;
                else
                    C_tmp = (max((h_UT-13)/10, 0)).^1.5;
                    Pr_LOS = (18/d_2Dout+exp(-d_2Dout/63)*(1-18/d_2Dout))*...
                            (1+C_tmp*(5/4)*((d_2Dout/100).^3)*exp(-d_2Dout/150));
                end
            elseif strcmp(obj.scenario.name,'3D-UMi')
                Pr_LOS = min(18./d_2Dout,1).*(1-exp(-d_2Dout./36))+exp(-d_2Dout./36);
            elseif strcmp(obj.scenario.name,'3D-RMa')
                Pr_LOS = min(exp(-(d_2Dout-10)./1000),1);
            elseif strcmp(obj.scenario.name,'3D-InH')
                d_2Din = obj.d_2D_in;
                switch obj.scenario.InH_case
                    case 'A'
                        Pr_LOS = min(max(exp(-(d_2Din-18)./27),0.5),1);
                    case {'B','open_office'}
                        Pr_LOS = min(max(exp(-(d_2Din-5)./70.8),0.54*exp(-(d_2Din-49)./211.7)),1);
                    case 'mixed_office'
                        Pr_LOS = min(max(exp(-(d_2Din-1.2)/4.7),0.32*exp(-(d_2Din-6.5)/32.6)),1);
                    otherwise
                        error(['Case "' obj.scenario.InH_case '" for 3D-InH is not valid!']);
                end
            end
        end
        
        % caculator the path loss distance depend
        function pathloss(obj)
            h_BS = obj.sector.h_BS;
            h_UT = obj.UE.h_UT;
            fcG = obj.fc/1e9;    % to GHz

            if strcmp(obj.scenario.name,'3D-UMa')
%                 W = obj.scenario.W;
%                 h = obj.scenario.h;
                if h_UT <= 13
                    C = 0;
                elseif (h_UT > 13) && (h_UT <= 23)
                    if obj.d_2D <= 18
                        g = 0;
                    else
                        g = 5/4*((obj.d_2D/100).^3)*exp(-obj.d_2D/150);
                    end
                    C = (((h_UT-13)/10).^1.5)*g;
                end
                if rand < 1/(1+C)
                    hE = 1;
                else
                    hE = 3*randi([4, obj.UE.n_fl],1);
%                     hE = randsrc(1,1,[12,15,18,21]);
                end
                h_BS_tmp = h_BS-hE;
                h_UT_tmp = h_UT-hE;
                d_BP_tmp = 4*h_BS_tmp*h_UT_tmp*obj.fc/3e8;
                if (obj.d_2D >= 10) && (obj.d_2D <= d_BP_tmp)
                    PL_LOS = 20*log10(obj.d_3D)+32.4+20*log10(fcG);
                elseif (obj.d_2D > d_BP_tmp) && (obj.d_2D < 5000)
                    PL_LOS = 40*log10(obj.d_3D)+32.4+20*log10(fcG)-10*log10(d_BP_tmp.^2+(h_BS-h_UT).^2);
                end
                if obj.bLOS
                    PL_db = PL_LOS;
                    sigma_SF_db = 4;
                else
                    % TR36.873
%                     PL_NLOS = 161.04-7.1*log10(W)+7.5*log10(h)-(24.37-3.7*(h/h_BS).^2)*log10(h_BS)+(43.42-3.1*log10(h_BS))*...
%                               (log10(obj.d_3D)-3)+20*log10(fcG)-(3.2*(log10(17.625)).^2-4.97)-0.6*(h_UT-1.5);
                    % TR38.901
                    PL_NLOS = 13.54 + 39.08*log10(obj.d_3D) + 20*log10(fcG)-0.6*(h_UT-1.5);
                    PL_db = max(PL_LOS, PL_NLOS);
                    sigma_SF_db = 6;
                    % optional
                    % PL_NLOS = 32.4+20*log10(fcG)+30*log10(obj.d_3D);
                    % sigma_SF_db = 7.8;
                end
                if obj.O2I
                    [P_tw,sigma_p] = obj.O2IpenetrationLoss(fcG,obj.UE.O2IPL);
                    sigma_SF_db = 7;
                    PL_db = PL_db + P_tw + 0.5*obj.d_2D_in + sigma_p*obj.UE.O2Isigma;
                end
            elseif strcmp(obj.scenario.name,'3D-UMi')
                h_BS_tmp = h_BS-1;
                h_UT_tmp = h_UT-1;
                d_BP_tmp = 4*h_BS_tmp*h_UT_tmp*obj.fc/3e8;
                if (obj.d_2D <= d_BP_tmp)% && (obj.d_2D >= 10)
                    % PL_LOS = 22*log10(obj.d_3D)+28+20*log10(fcG); TR36.873
                    PL_LOS = 32.4 + 21*log10(obj.d_3D) + 20*log10(fcG); % TR38.901
                elseif  (obj.d_2D > d_BP_tmp) && (obj.d_2D < 5000)
                    % PL_LOS = 40*log10(obj.d_3D)+28+20*log10(fcG)-9*log10(d_BP_tmp.^2+(h_BS-h_UT).^2); TR36.873
                    PL_LOS = 32.4 + 40*log10(obj.d_3D) + 20*log10(fcG) - 9.5*log10(d_BP_tmp.^2+(h_BS-h_UT).^2); % TR38.901
                end
                if obj.bLOS
                    PL_db = PL_LOS;
                    % sigma_SF_db = 3; TR36.873
                    sigma_SF_db = 4; % TR38.901
                else
                    % PL_NLOS = 36.7*log10(obj.d_3D)+22.7+26*log10(fcG)-0.3*(h_UT-1.5); TR36.873
                    % sigma_SF_db = 4; TR36.873
                    PL_NLOS = 35.3*log10(obj.d_3D)+22.4+21.3*log10(fcG)-0.3*(h_UT-1.5); % TR38.901
                    PL_db = max(PL_LOS, PL_NLOS);
                    sigma_SF_db = 7.82; % TR38.901
                    
                    % optional
                    % PL_NLOS = 32.4+20*log10(fcG)+31.9*log10(obj.d_3D);
                    % sigma_SF_db = 8.2;
                end
                if obj.O2I
                    [P_tw,sigma_p] = obj.O2IpenetrationLoss(fcG,obj.UE.O2IPL);
                    sigma_SF_db = 7;
                    PL_db = PL_db + P_tw + 0.5*obj.d_2D_in + sigma_p*obj.UE.O2Isigma;
                end
            elseif strcmp(obj.scenario.name,'3D-RMa')
                W = obj.scenario.W;
                h = obj.scenario.h;
                d_BP = (2*pi*h_BS*h_UT*obj.fc)/3e8;
                PL1 = @(d_3D)20*log10(40*pi*d_3D*fcG/3) + min(0.03*(h^1.72),10)*log10(d_3D)...
                      -min(0.044*(h^1.72), 14.77) + 0.002*log10(h)*d_3D;
                if (obj.d_2D > 10) && (obj.d_2D <= d_BP)
                    PL_LOS = PL1(obj.d_3D);
                    sigma_SF_db = 4;
                elseif  (obj.d_2D > d_BP) && (obj.d_2D <= 10000)
                    PL_LOS = PL1(d_BP) + 40*log10(obj.d_3D/d_BP);
                    sigma_SF_db = 6;
                end
                if obj.bLOS
                    PL_db = PL_LOS;
                else
                    PL_db = 161.04-7.1*log10(W)+7.5*log10(h)-(24.37-3.7*(h/h_BS)^2)*log10(h_BS)+ ...
                            (43.42-3.1*log10(h_BS))*(log10(obj.d_3D)-3)+20*log10(fcG)-(3.2*(log10(11.75*h_UT))^2-4.97);
                    sigma_SF_db = 8;
                end
                if obj.O2I
                    [P_tw,sigma_p] = obj.O2IpenetrationLoss(fcG,obj.UE.O2IPL);
                    sigma_SF_db = 8;
                    PL_db = PL_db + P_tw + 0.5*obj.d_2D_in + sigma_p*obj.UE.O2Isigma;
                end
            elseif strcmp(obj.scenario.name,'3D-InH')
                % TR36.873
                % PL_db = 16.9*log10(obj.d_3D) + 32.8 + 20*log10(fcG);
                % sigma_SF_db = 3;
                % TR38.901
                PL_db = 32.4 + 17.3*log10(obj.d_3D) + 20*log10(fcG);
                sigma_SF_db = 3;
                if ~obj.bLOS
                    % TR36.873
                    % PL_db = 43.3*log10(obj.d_3D) + 11.5 + 20*log10(fcG);
                    % sigma_SF_db = 4;
                    % TR38.901
                    PL_NLOS = 17.3 + 38.3*log10(obj.d_3D) + 24.9*log10(fcG);
                    PL_db = max(PL_db,PL_NLOS);
                    sigma_SF_db = 8.03;
                    % optional
                    % PL_NLOS = 32.4+20*log10(fcG)+31.9*log10(obj.d_3D);
                    % sigma_SF_db = 8.29;
                end
                
%             elseif strcmp(obj.scenario.name,'InF')
%                 PL_db = 31.84+21.5*log10(obj.d_3D)+19*log10(fcG);
%                 sigma_SF_db = 4;
%                 if ~obj.bLOS
%                     switch obj.scenario.InF_case
%                         case 'SL'
%                             PL_db = max(PL_db,33+25.5*log10(obj.d_3D)+20*log10(fcG));
%                             sigma_SF_db = 5.7;
%                         case 'DL'
%                             PL_db = max(PL_db,18.6+35.7*log10(obj.d_3D)+20*log10(fcG));
%                             sigma_SF_db = 7.2;
%                         case 'SH'
%                             PL_db = max(PL_db,32.4+23*log10(obj.d_3D)+20*log10(fcG));
%                             sigma_SF_db = 5.9;
%                         case 'DH'
%                             PL_db = max(PL_db,33.63+21.9*log10(obj.d_3D)+20*log10(fcG));
%                     end
%                 end
            end          
            
            obj.PL = PL_db;
            obj.sigma_SF = sigma_SF_db;
        end
        
        function [P_tw,sigma_p] = O2IpenetrationLoss(obj,fc,LowHigh)
            if fc<6
                P_tw = 20; sigma_p = 0;
            else
                if strcmp(obj.scenario.name,'3D-RMa')
                    LowHigh = 'low';
                elseif strcmp(obj.scenario.name,'InF')
                    LowHigh = 'high';
                end
                MaterialLoss = tools.get_MaterialLoss(fc);
                if strcmp(LowHigh,'low')
                    P_tw = 5-10*log10(0.3*10.^(-MaterialLoss(1,:)/10) + 0.7*10.^(-MaterialLoss(3,:)/10));
                    sigma_p = 4.4;
                elseif strcmp(LowHigh,'high')
                    P_tw = 5-10*log10(0.7*10.^(-MaterialLoss(2,:)/10) + 0.3*10.^(-MaterialLoss(3,:)/10));
                    sigma_p = 6.5;
                end
            end
        end
        
        % caculator the large scale parameters
        function large_scale_para(obj, tab)
            fc = obj.fc/1e9; %#ok<*PROPLC>
            if strcmp(obj.scenario.name,'3D-UMa')
                fc(fc<=6) = 6;
                if obj.O2I
                    if obj.bLOS
                        obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)-0.01*(obj.UE.h_UT-1.5)+0.75);
                        obj.sigma_lgZSD = 0.4;
                        obj.mu_offset_ZOD = 0;
                    else
                        obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)-0.01*(obj.UE.h_UT-1.5)+0.9);
                        obj.sigma_lgZSD = 0.49;
                        % obj.mu_offset_ZOD = -10^(-0.62*log10(max(10, obj.d_2D_out))+1.93-0.07*(obj.UE.h_UT-1.5));
                        obj.mu_offset_ZOD = (7.66*log10(fc)-5.96) -10^((0.208*log10(fc)-0.782)*log10(max((25), obj.d_2D_out))+(-0.13*log10(fc)+2.03));
                    end
                    temp1 = obj.sector.corr_sos{3}.randn([obj.UE.pos,obj.UE.h_UT]'); % obj.UE.n_fl
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.UMa_C_sqrt_O2I*ksi;
                    % SF/DS/ASD/ASA/ZSD/ZSA
                    std_error = [obj.sigma_SF, 0.32, 0.42, 0.16, obj.sigma_lgZSD, 0.43]';
                    result = std_error.*ksi;
                    obj.SF = result(1);
                    obj.DS = 10^(result(2)-6.63);
                    obj.ASD = 10^(result(3)+1.25);
                    obj.ASA = 10^(result(4)+1.76);
                    obj.ZSD = 10^(result(5)+obj.mu_lgZSD);
                    obj.ZSA = 10^(result(6)+1.01);
                    obj.r_tau = 2.2;
                    obj.mu_XPR = 9;
                    obj.sigma_XPR = 5;  % 5
                    obj.N = 12;
                    obj.M = 20;
                    obj.zeta = 4;
                    obj.c_DS  = 11;
                    obj.c_ASA = 20;
                    obj.c_ZSA = 6;
                    obj.c_ASD = 5;
                elseif obj.bLOS % 0 1
                    obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.UE.h_UT-1.5)+0.75);
                    obj.sigma_lgZSD = 0.4;
                    obj.mu_offset_ZOD = 0;
                    temp1 = obj.sector.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.UMa_C_sqrt_LOS*ksi;
                    % SF/K/DS/ASD/ASA/ZSD/ZSA
                    std_error = [obj.sigma_SF, 3.5, 0.66, 0.28, 0.20, obj.sigma_lgZSD, 0.16]';
                    result = std_error.*ksi;
                    obj.SF = result(1);
                    obj.K = result(2)+9;
                    obj.DS = 10^(result(3)+(-6.955-0.0963*log10(fc)));
                    obj.ASD = 10^(result(4)+(1.06+0.1114*log10(fc)));
                    obj.ASA = 10^(result(5)+1.81);
                    obj.ZSD = 10^(result(6)+obj.mu_lgZSD);
                    obj.ZSA = 10^(result(7)+0.95);
                    obj.r_tau = 2.5;
                    obj.mu_XPR = 8;
                    obj.sigma_XPR = 4;
                    obj.N = 12;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_DS  = max(0.25,6.5622-3.4084*log10(fc));
                    obj.c_ASA = 11;
                    obj.c_ZSA = 7;
                    obj.c_ASD = 5;
                else % 1 0
%                     weight_bs = 1; weight_ue = 0;
                    obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.UE.h_UT-1.5)+0.9);
                    obj.sigma_lgZSD = 0.49;
                    % obj.mu_offset_ZOD = -10^(-0.62*log10(max(10, obj.d_2D))+1.93-0.07*(obj.UE.h_UT-1.5));
                    obj.mu_offset_ZOD = (7.66*log10(fc)-5.96) -10^((0.208*log10(fc)-0.782)*log10(max((25), obj.d_2D))+(-0.13*log10(fc)+2.03));
                    temp1 = obj.sector.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.UMa_C_sqrt_NLOS*ksi;
                    % SF/DS/ASD/ASA/ZSD/ZSA
                    std_error = [obj.sigma_SF, 0.39, 0.28, 0.11, obj.sigma_lgZSD, 0.16]';
                    result = std_error.*ksi;
                    obj.SF = result(1);   % dB
                    obj.DS = 10^(result(2)+(-6.28-0.204*log10(fc)));
                    obj.ASD = 10^(result(3)+(1.5-0.1144*log10(fc)));
                    obj.ASA = 10^(result(4)+(2.08-0.27*log10(fc)));
                    obj.ZSD = 10^(result(5)+obj.mu_lgZSD);
                    obj.ZSA = 10^(result(6)+(-0.3236*log10(fc)+1.512));
                    obj.r_tau = 2.3;
                    obj.mu_XPR = 7;
                    obj.sigma_XPR = 3;
                    obj.N = 20;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_DS  = max(0.25,6.5622-3.4084*log10(fc));
                    obj.c_ASA = 15;
                    obj.c_ZSA = 7;
                    obj.c_ASD = 2;
                end
            elseif strcmp(obj.scenario.name,'3D-UMi')
                fc(fc<=2) = 2;
                if obj.O2I % 
                    if obj.bLOS
                        % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*abs(obj.UE.h_UT-obj.sector.attached_BS.h_BS)+0.75);
                        obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D/1000)+0.01*abs(obj.UE.h_UT-obj.sector.h_BS)+0.83);
                        obj.sigma_lgZSD = 0.35; % 0.4
                        obj.mu_offset_ZOD = 0;
                    else
                        % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*max(obj.UE.h_UT-obj.sector.attached_BS.h_BS, 0)+0.9);
                        obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D/1000)+0.01*max(obj.UE.h_UT-obj.sector.h_BS, 0)+0.2);
                        obj.sigma_lgZSD = 0.35; %0.6
                        % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D_out))+1.6);
                        obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D))+3.3);
                    end
                    temp1 = obj.sector.corr_sos{3}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.UMi_C_sqrt_O2I*ksi;
                    std_error = [obj.sigma_SF, 0.32, 0.42, 0.16, obj.sigma_lgZSD, 0.43]';
                    result = std_error.*ksi;
                    obj.SF = result(1);
                    obj.DS = 10^(result(2)-6.62);
                    obj.ASD = 10^(result(3)+1.25);
                    obj.ASA = 10^(result(4)+1.76);
                    obj.ZSD = 10^(result(5)+obj.mu_lgZSD);
                    obj.ZSA = 10^(result(6)+1.01);
                    obj.r_tau = 2.2;
                    obj.mu_XPR = 9;
                    obj.sigma_XPR = 5;  % 5
                    obj.N = 12;
                    obj.M = 20;
                    obj.zeta = 4;
                    obj.c_DS  = 11;
                    obj.c_ASA = 20;
                    obj.c_ZSA = 6;
                    obj.c_ASD = 5;
                elseif obj.bLOS % 0 1
                    % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*abs(obj.UE.h_UT-obj.sector.attached_BS.h_BS)+0.75);
                    obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D/1000)+0.01*abs(obj.UE.h_UT-obj.sector.h_BS)+0.83);
                    obj.sigma_lgZSD = 0.35; % 0.4
                    obj.mu_offset_ZOD = 0;
                    temp1 = obj.sector.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.UMi_C_sqrt_LOS*ksi;
                    % SF/K/DS/ASD/ASA/ZSD/ZSA
                    % std_error = [obj.sigma_SF, 5, 0.4, 0.43, 0.19, obj.sigma_lgZSD, 0.16]';
                    std_error = [obj.sigma_SF, 5, 0.38, 0.41, (0.014*log10(1+fc)+0.28), obj.sigma_lgZSD, (-0.04*log10(1+fc)+0.34)]';
                    result = std_error.*ksi;
                    obj.SF = result(1);
                    obj.K = result(2)+9;
                    obj.DS = 10^(result(3)+(-0.24*log10(1+fc)-7.14));
                    obj.ASD = 10^(result(4)+(-0.05*log10(1+fc)+1.21));
                    obj.ASA = 10^(result(5)+(-0.08*log10(1+fc)+1.73));
                    obj.ZSD = 10^(result(6)+obj.mu_lgZSD);
                    obj.ZSA = 10^(result(7)+(-0.1*log10(1+fc)+0.73));
                    obj.r_tau = 3;
                    obj.mu_XPR = 9;
                    obj.sigma_XPR = 3;
                    obj.N = 12;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_DS  = 5;
                    obj.c_ASA = 17;
                    obj.c_ZSA = 7;
                    obj.c_ASD = 3;
                else % NLOS
                    % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*max(obj.UE.h_UT-obj.sector.attached_BS.h_BS, 0)+0.9);
                    obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D/1000)+0.01*max(obj.UE.h_UT-obj.sector.h_BS, 0)+0.2);
                    obj.sigma_lgZSD = 0.35;
                    % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D))+1.6);
                    obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D))+3.3);
                    temp1 = obj.sector.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.UMi_C_sqrt_NLOS*ksi;
                    % SF/DS/ASD/ASA/ZSD/ZSA
                    std_error = [obj.sigma_SF, (0.16*log10(1+fc)+0.28), (0.11*log10(1+fc)+0.33), (0.05*log10(1+fc)+0.3), obj.sigma_lgZSD, (-0.07*log10(1+fc)+0.41)]';
                    result = std_error.*ksi;
                    obj.SF = result(1);
                    obj.DS = 10^(result(2)+(-0.24*log10(1+fc)-6.83));
                    obj.ASD = 10^(result(3)+(-0.23*log10(1+fc)+1.53));
                    obj.ASA = 10^(result(4)+(-0.08*log10(1+fc)+1.81));
                    obj.ZSD = 10^(result(5)+obj.mu_lgZSD);
                    obj.ZSA = 10^(result(6)+(-0.04*log10(1+fc)+0.92));
                    obj.r_tau = 2.1;
                    obj.mu_XPR = 8;
                    obj.sigma_XPR = 3;
                    obj.N = 19;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_DS  = 11;
                    obj.c_ASA = 22;
                    obj.c_ZSA = 7;
                    obj.c_ASD = 10;
                end
            elseif strcmp(obj.scenario.name,'3D-RMa') % without calibration
                if obj.O2I
                    if obj.bLOS
                        obj.mu_lgZSD = 0.3;
                        obj.sigma_lgZSD = 0.4;
                        obj.mu_offset_ZOD = 0;
                        temp1 = obj.sector.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
                        ksi = reshape(temp1,[numel(temp1), 1]);
                        ksi = tab.RMa_C_sqrt_LOS*ksi;
                        obj.SF = obj.sigma_SF*ksi(1);
                        obj.K = 4*ksi(2)+7;
                        obj.DS = 10^(0.55*ksi(3)-7.49);
                        obj.ASD = 10^(0.38*ksi(4)+0.90);
                        obj.ASA = 10^(0.24*ksi(5)+1.52);
                        obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
                        obj.ZSA = 10^(0.16*ksi(7)+0.60);
                        obj.r_tau = 3.8;
                        obj.mu_XPR = 12;
                        obj.sigma_XPR = 4;
                        obj.N = 11;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_ASA = 3;
                        obj.c_ZSA = 3;
                        obj.c_ASD = 2;
                    else
                        obj.mu_lgZSD = 0.3;
                        obj.sigma_lgZSD = 0.49;
                        obj.mu_offset_ZOD = atand((35-5)/obj.d_2D)-atand((35-1.5)/obj.d_2D);
                        temp1 = obj.sector.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
                        ksi = reshape(temp1,[numel(temp1), 1]);
                        ksi = tab.RMa_C_sqrt_NLOS*ksi;
                        obj.SF = obj.sigma_SF*ksi(1);
                        obj.DS = 10^(0.48*ksi(2)-7.43);
                        obj.ASD = 10^(0.45*ksi(3)+0.95);
                        obj.ASA = 10^(0.13*ksi(4)+1.52);
                        obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
                        obj.ZSA = 10^(0.16*ksi(6)+0.88);
                        obj.r_tau = 1.7;
                        obj.mu_XPR = 7;
                        obj.sigma_XPR = 3;
                        obj.N = 10;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_ASA = 3;
                        obj.c_ZSA = 3;
                        obj.c_ASD = 2;
                    end
                elseif obj.bLOS
                    obj.mu_lgZSD = 0.3;
                    obj.sigma_lgZSD = 0.4;
                    obj.mu_offset_ZOD = 0;
                    temp1 = obj.sector.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.RMa_C_sqrt_LOS*ksi;
                    obj.SF = obj.sigma_SF*ksi(1);
                    obj.K = 4*ksi(2)+7;
                    obj.DS = 10^(0.55*ksi(3)-7.49);
                    obj.ASD = 10^(0.38*ksi(4)+0.90);
                    obj.ASA = 10^(0.24*ksi(5)+1.52);
                    obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
                    obj.ZSA = 10^(0.16*ksi(7)+0.6);
                    obj.r_tau = 3.8;
                    obj.mu_XPR = 12;
                    obj.sigma_XPR = 4;
                    obj.N = 11;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_ASA = 3;
                    obj.c_ZSA = 3;
                    obj.c_ASD = 2;
                else
                    obj.mu_lgZSD = 0.3;
                    obj.sigma_lgZSD = 0.49;
                    obj.mu_offset_ZOD = atand((35-5)/obj.d_2D)-atand((35-1.5)/obj.d_2D);
                    temp1 = obj.sector.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.RMa_C_sqrt_NLOS*ksi;
                    obj.SF = obj.sigma_SF*ksi(1);
                    obj.DS = 10^(0.48*ksi(2)-7.43);
                    obj.ASD = 10^(0.45*ksi(3)+0.95);
                    obj.ASA = 10^(0.13*ksi(4)+1.52);
                    obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
                    obj.ZSA = 10^(0.16*ksi(6)+0.88);
                    obj.r_tau = 1.7;
                    obj.mu_XPR = 7;
                    obj.sigma_XPR = 3;
                    obj.N = 10;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_ASA = 3;
                    obj.c_ZSA = 3;
                    obj.c_ASD = 2;
                end
            elseif  strcmp(obj.scenario.name,'3D-InH')
                fc(fc<=6) = 6;
                if obj.bLOS
                    % obj.mu_lgZSD = 1.02;
                    % obj.sigma_lgZSD = 0.41;
                    obj.mu_lgZSD = -1.43*log10(1+fc)+2.25;
                    obj.sigma_lgZSD = 0.13*log10(1+fc)+0.15;
                    obj.mu_offset_ZOD = 0;
                    temp1 = obj.sector.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.InH_C_sqrt_LOS*ksi;
                    obj.SF = obj.sigma_SF*ksi(1);
                    obj.K = (-0.58*log10(1+fc)+6.19)*ksi(2)+ (0.84*log10(1+fc)+2.12);
                    obj.DS = 10^((-0.16*log10(1+fc)+0.5)*ksi(3)+ (-0.01*log10(1+fc)-7.79));
                    obj.ASD = 10^(0.18*ksi(4)+1.60);
                    obj.ASA = 10^((0.12*log10(1+fc))*ksi(5)+ (-0.19*log10(1+fc)+1.86));
                    obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
                    obj.ZSA = 10^((-0.04*log10(1+fc)+0.17)*ksi(7)+ (-0.26*log10(1+fc)+1.21));
                    obj.r_tau = 2.15;
                    obj.mu_XPR = 15;
                    obj.sigma_XPR = 3; % 3
                    obj.N = 8;
                    obj.M = 20;
                    obj.zeta = 6;
                    obj.c_ASA = -6.2*log10(1+fc)+16.72;
                    obj.c_ZSA = -3.85*log10(1+fc)+10.28;
                    obj.c_ASD = 7;
                else
                    obj.mu_lgZSD = 1.37;
                    obj.sigma_lgZSD = 0.38;
                    obj.mu_offset_ZOD = 0;
                    temp1 = obj.sector.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
                    ksi = reshape(temp1,[numel(temp1), 1]);
                    ksi = tab.InH_C_sqrt_NLOS*ksi;
                    obj.SF = obj.sigma_SF*ksi(1);
                    obj.DS = 10^((0.1*log10(1+fc)+0.11)*ksi(2)+ (-0.28*log10(1+fc)-7.29));
                    obj.ASD = 10^(0.17*ksi(3)+1.49);
                    obj.ASA = 10^((0.12*log10(1+fc))*ksi(4)+ (-0.11*log10(1+fc)+1.8));
                    obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
                    obj.ZSA = 10^((-0.09*log10(1+fc)+0.24)*ksi(6)+ (-0.15*log10(1+fc)+1.04));
                    obj.r_tau = 1.84;
                    obj.mu_XPR = 12;
                    obj.sigma_XPR = 7; % 3
                    obj.N = 10;
                    obj.M = 20;
                    obj.zeta = 3;
                    obj.c_ASA = -13*log10(1+fc)+30.53;
                    obj.c_ZSA = -3.72*log10(1+fc)+10.25;
                    obj.c_ASD = 3;
                end
%             elseif  strcmp(obj.scenario.name,'InF')
%                 % V = hall volume in m^3, S = total surface area of hall in m^2 (walls+floor+ceiling)
%                 if obj.bLOS
%                     obj.mu_lgZSD = 1.35;
%                     obj.sigma_lgZSD = 0.35;
%                     obj.mu_offset_ZOD = 0;
%                     temp1 = obj.sector.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
%                     ksi = reshape(temp1,[numel(temp1), 1]);
%                     ksi = tab.InH_C_sqrt_LOS*ksi;
%                     obj.SF = obj.sigma_SF*ksi(1);
%                     obj.K = 8*ksi(2)+7;
%                     obj.DS = 10^(0.15*ksi(3)+ (log10(26*(obj.scenario.V/obj.scenario.S)+14)-9.35));
%                     obj.ASD = 10^(0.25*ksi(4)+1.56);
%                     obj.ASA = 10^((0.12*log10(1+fc)+0.2)*ksi(5)+ (-0.18*log10(1+fc)+1.78));
%                     obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
%                     obj.ZSA = 10^(0.35*ksi(7)+ (-0.2*log10(1+fc)+1.5));
%                     obj.r_tau = 2.7;
%                     obj.mu_XPR = 12;
%                     obj.sigma_XPR = 6;
%                     obj.N = 25;
%                     obj.M = 20;
%                     obj.zeta = 4;
%                     obj.c_ASA = 8;
%                     obj.c_ZSA = 9;
%                     obj.c_ASD = 5;
%                 else
%                     obj.mu_lgZSD = 1.2;
%                     obj.sigma_lgZSD = 0.55;
%                     obj.mu_offset_ZOD = 0;
%                     temp1 = obj.sector.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
%                     ksi = reshape(temp1,[numel(temp1), 1]);
%                     ksi = tab.InH_C_sqrt_NLOS*ksi;
%                     obj.SF = obj.sigma_SF*ksi(1);
%                     obj.DS = 10^(0.19*ksi(2)+ (log10(30*(obj.scenario.V/obj.scenario.S)+32)-9.44));
%                     obj.ASD = 10^(0.2*ksi(3)+1.57);
%                     obj.ASA = 10^(0.3*ksi(4)+1.72);
%                     obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
%                     obj.ZSA = 10^(0.45*ksi(6)+ (-0.13*log10(1+fc)+1.45));
%                     obj.r_tau = 3;
%                     obj.mu_XPR = 11;
%                     obj.sigma_XPR = 6;
%                     obj.N = 25;
%                     obj.M = 20;
%                     obj.zeta = 3;
%                     obj.c_ASA = 8;
%                     obj.c_ZSA = 9;
%                     obj.c_ASD = 5;
%                 end
            end
        end
        
        % caculator the AOA for every ray
        function AOA_calc(obj)
%             if strcmp(obj.scenario.name,'3D-InH')
%                 C_phi_NLOS_TAB = [1.434, 1,1,1, 1.501];
%                 C_phi_NLOS = C_phi_NLOS_TAB(obj.N - 14);
%                 if obj.bLOS && ~obj.O2I
%                     C_phi = C_phi_NLOS*(0.9275-0.0439*obj.K-0.0071*obj.K^2+0.0002*obj.K^3);
%                 else 
%                     C_phi = C_phi_NLOS;
%                 end
%                 phi_n_AOA_tmp = -obj.ASA*log(obj.Pn_LOS/max(obj.Pn_LOS))/C_phi;
%             else  % for 3D-UMa, 3D-RMa, 3D-UMi
%                 %                                    4       5                8           10     11     12           14     15     16                    19     20
%                 C_phi_NLOS_TAB = [0.5, 0.58, 0.69, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
% %                 C_phi_NLOS_TAB = [0.779, 0.779, 0.779, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
%                 C_phi_NLOS = C_phi_NLOS_TAB(obj.N);
%                 if obj.bLOS && ~obj.O2I
%                     C_phi = C_phi_NLOS*(1.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3);
%                 else
%                     C_phi = C_phi_NLOS;
%                 end
%                 phi_n_AOA_tmp = 2*(obj.ASA/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/C_phi;
%             end
            n_clusters     = [4 5 8 10 11 12 14 15 16 19 20 25];
            C_phi_NLOS_TAB = [0.779, 0.860, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.358];
            C_phi_NLOS     = C_phi_NLOS_TAB(n_clusters == obj.N);
            C_phi          = sum([C_phi_NLOS ,C_phi_NLOS*(0.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3)]); % for NLOS and O2I, K = [].
            phi_n_AOA_tmp  = 2*(obj.ASA/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/C_phi;
            
            if ~isempty(obj.sector.SpatialConsistency)
                Xn = round(obj.sector.SpatialConsistency.AOAsign{obj.spatialConsist_indx}.rand(obj.UE.initPos))*2-1;
                Yn = (obj.ASA/7)*obj.sector.SpatialConsistency.AOAoffset{obj.spatialConsist_indx}.randn([obj.UE.pos,obj.UE.h_UT]');
                Xn = Xn(obj.tau_order);
                Yn = Yn(obj.tau_order);
                Xn = Xn(obj.keep).';
                Yn = Yn(obj.keep).';
            else
                Xn = randsrc(1,obj.N_new,[-1, 1]);
                Yn = (obj.ASA/7)*randn(1,obj.N_new);
            end
            phi_n_AOA = Xn.*phi_n_AOA_tmp+Yn+obj.phi_LOS_AOA;
            if obj.bLOS % && ~obj.O2I
                phi_n_AOA = (Xn.*phi_n_AOA_tmp+Yn)-(Xn(1)*phi_n_AOA_tmp(1)+Yn(1)-obj.phi_LOS_AOA);
            end
            alpha_m = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                       0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            phi_n_m_AOA_ = repmat(phi_n_AOA.',1,obj.M)+obj.c_ASA(ones(obj.N_new,1))*alpha_m;
            phi_n_m_AOA_ = mod(phi_n_m_AOA_,360);
            phi_n_m_AOA_ = phi_n_m_AOA_-360*floor(phi_n_m_AOA_/180);
            obj.phi_n_m_AOA = phi_n_m_AOA_;
        end
        
        % caculator the AOD for every ray
        function AOD_calc(obj)
%             if strcmp(obj.scenario.name,'3D-InH')
%                 C_phi_NLOS_TAB = [1.434, 1,1,1, 1.501];
%                 C_phi_NLOS = C_phi_NLOS_TAB(obj.N - 14);
%                 if obj.bLOS && ~obj.O2I 
%                     C_phi = C_phi_NLOS*(0.9275-0.0439*obj.K-0.0071*obj.K^2+0.0002*obj.K^3);
%                 else 
%                     C_phi = C_phi_NLOS;
%                 end
%                 phi_n_AOD_tmp = -obj.ASD*log(obj.Pn_LOS/max(obj.Pn_LOS))/C_phi;
%             else  % for 3D-UMa, 3D-RMa, 3D-UMi
%                 %                                    4       5                8           10     11     12           14     15     16                    19     20
%                 C_phi_NLOS_TAB = [0.5, 0.58, 0.69, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
% %                 C_phi_NLOS_TAB = [0.779, 0.779, 0.779, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
%                 C_phi_NLOS = C_phi_NLOS_TAB(obj.N);
%                 if obj.bLOS && ~obj.O2I
%                     C_phi = C_phi_NLOS*(1.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3);
%                 else
%                     C_phi = C_phi_NLOS;
%                 end
%                 phi_n_AOD_tmp = 2*(obj.ASD/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/C_phi;
%             end
            n_clusters     = [4 5 8 10 11 12 14 15 16 19 20 25];
            C_phi_NLOS_TAB = [0.779, 0.860, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.358];
            C_phi_NLOS     = C_phi_NLOS_TAB(n_clusters == obj.N);
            C_phi          = sum([C_phi_NLOS ,C_phi_NLOS*(0.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3)]); % for NLOS and O2I, K = [].
            phi_n_AOD_tmp  = 2*(obj.ASD/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/C_phi;
%             phi_n_AOD_tmp  = -obj.ASD*log(obj.Pn_LOS/max(obj.Pn_LOS))/C_phi;
            
            if ~isempty(obj.sector.SpatialConsistency)
                Xn = round(obj.sector.SpatialConsistency.AODsign{obj.spatialConsist_indx}.rand([obj.sector.initalPOS, obj.sector.h_BS]'))*2-1;
                Yn = (obj.ASD/7)*obj.sector.SpatialConsistency.AODoffset{obj.spatialConsist_indx}.randn([obj.sector.Position, obj.sector.h_BS]');
                Xn = Xn(obj.tau_order);
                Yn = Yn(obj.tau_order);
                Xn = Xn(obj.keep).';
                Yn = Yn(obj.keep).';
            else
                Xn = randsrc(1,obj.N_new,[-1, 1]);  
                Yn = (obj.ASD/7)*randn(1,obj.N_new);
            end
            phi_n_AOD = Xn.*phi_n_AOD_tmp+Yn+obj.phi_LOS_AOD;
            if obj.bLOS % && ~obj.O2I
                phi_n_AOD = (Xn.*phi_n_AOD_tmp+Yn)-(Xn(1)*phi_n_AOD_tmp(1)+Yn(1)-obj.phi_LOS_AOD);
            end
            alpha_m = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                       0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            phi_n_m_AOD_ = repmat(phi_n_AOD.',1,obj.M)+obj.c_ASD(ones(obj.N_new,1))*alpha_m;
            phi_n_m_AOD_ = mod(phi_n_m_AOD_,360);
            phi_n_m_AOD_ = phi_n_m_AOD_-360*floor(phi_n_m_AOD_/180);
            obj.phi_n_m_AOD = phi_n_m_AOD_;
        end
        
        % caculator the ZOA for every ray
        function ZOA_calc(obj)
%             %                                                                         10      11    12                    15                            19     20
%             C_phi_NLOS_TAB = [0.582,0.642,0.698,0.750,0.798,0.842,0.882,0.92,0.954, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
% %             C_phi_NLOS_TAB = [0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
%             C_phi_NLOS = C_phi_NLOS_TAB(obj.N);
            n_clusters       = [8 10 11 12 15 19 20 25];
            C_theta_NLOS_TAB = [0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178, 1.282];
            C_theta_NLOS     = C_theta_NLOS_TAB(n_clusters == obj.N);
            C_theta          = sum([C_theta_NLOS ,C_theta_NLOS*(0.3086+0.0339*obj.K-0.0077*obj.K^2+0.0002*obj.K^3)]); % for NLOS and O2I, K = [].
            theta_n_ZOA_tmp  = -obj.ZSA*log(obj.Pn_LOS/max(obj.Pn_LOS))/C_theta;
            
            if ~isempty(obj.sector.SpatialConsistency)
                Xn = round(obj.sector.SpatialConsistency.ZOAsign{obj.spatialConsist_indx}.rand([obj.sector.initalPOS, obj.sector.h_BS]'))*2-1;
                Yn = (obj.ZSA/7)*obj.sector.SpatialConsistency.ZOAoffset{obj.spatialConsist_indx}.randn([obj.UE.pos,obj.UE.h_UT]');
                Xn = Xn(obj.tau_order);
                Yn = Yn(obj.tau_order);
                Xn = Xn(obj.keep).';
                Yn = Yn(obj.keep).';
            else
                Xn  = randsrc(1,obj.N_new,[-1, 1]);  
                Yn  = (obj.ZSA/7)*randn(1,obj.N_new);
            end
            ZOA = obj.theta_LOS_ZOA;
            ZOA(obj.O2I) = 90;
            theta_n_ZOA = Xn.*theta_n_ZOA_tmp+Yn+ZOA;
            if obj.bLOS  && ~obj.O2I
                theta_n_ZOA = (Xn.*theta_n_ZOA_tmp+Yn)-(Xn(1)*theta_n_ZOA_tmp(1)+Yn(1)-obj.theta_LOS_ZOA);
            end
            alpha_m = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                       0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            theta_n_m_ZOA_ = repmat(theta_n_ZOA.',1,obj.M) + obj.c_ZSA(ones(obj.N_new,1))*alpha_m;
            theta_n_m_ZOA_t = mod((theta_n_m_ZOA_+360),360);  % 360
            index = ((theta_n_m_ZOA_t <= 360)&(theta_n_m_ZOA_t >= 180));
            theta_n_m_ZOA_ = (~index).*theta_n_m_ZOA_t + index.*(360 - theta_n_m_ZOA_t);
            obj.theta_n_m_ZOA = theta_n_m_ZOA_;
        end
        
        % caculator the ZOD for every ray
        function ZOD_calc(obj)
%             %                                                                         10      11    12                    15                            19     20
%             C_phi_NLOS_TAB = [0.582,0.642,0.698,0.750,0.798,0.842,0.882,0.92,0.954, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
% %             C_phi_NLOS_TAB = [0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
%             C_phi_NLOS = C_phi_NLOS_TAB(obj.N);
            n_clusters       = [8 10 11 12 15 19 20 25];
            C_theta_NLOS_TAB = [0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178, 1.282];
            C_theta_NLOS     = C_theta_NLOS_TAB(n_clusters == obj.N);
            C_theta          = sum([C_theta_NLOS ,C_theta_NLOS*(0.3086+0.0339*obj.K-0.0077*obj.K^2+0.0002*obj.K^3)]); % for NLOS and O2I, K = [].
            theta_n_ZOD_tmp  = -obj.ZSD*log(obj.Pn_LOS/max(obj.Pn_LOS))/C_theta;
            
            if ~isempty(obj.sector.SpatialConsistency)
                Xn = round(obj.sector.SpatialConsistency.ZODsign{obj.spatialConsist_indx}.rand(obj.UE.initPos))*2-1;
                Yn = (obj.ZSD/7)*obj.sector.SpatialConsistency.ZODoffset{obj.spatialConsist_indx}.randn([obj.sector.Position, obj.sector.h_BS]');
                Xn = Xn(obj.tau_order);
                Yn = Yn(obj.tau_order);
                Xn = Xn(obj.keep).';
                Yn = Yn(obj.keep).';
            else
                Xn = randsrc(1,obj.N_new,[-1, 1]);  
                Yn = (obj.ZSD/7)*randn(1,obj.N_new);
            end
            theta_n_ZOD = Xn.*theta_n_ZOD_tmp+Yn+obj.theta_LOS_ZOD+obj.mu_offset_ZOD;
            if obj.bLOS % && ~obj.O2I
                theta_n_ZOD = (Xn.*theta_n_ZOD_tmp+Yn)-(Xn(1)*theta_n_ZOD_tmp(1)+Yn(1)-obj.theta_LOS_ZOD);
            end
            alpha_m = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                       0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            c_ZSD = (3/8)*(10^obj.mu_lgZSD);
            theta_n_m_ZOD_  = repmat(theta_n_ZOD.',1,obj.M) + c_ZSD(ones(obj.N_new,1))*alpha_m;
            theta_n_m_ZOD_t = mod((theta_n_m_ZOD_+360),360);   % 360
            index = ((theta_n_m_ZOD_t <= 360)&(theta_n_m_ZOD_t >= 180));
            theta_n_m_ZOD_ = (~index).*theta_n_m_ZOD_t + index.*(360 - theta_n_m_ZOD_t);
            obj.theta_n_m_ZOD = theta_n_m_ZOD_;
        end
        
        function RandomCouplingRays(obj)
            if ~isempty(obj.sector.SpatialConsistency)
                rn1 = obj.sector.SpatialConsistency.randomCoupling1{obj.spatialConsist_indx}.rand(obj.UE.initPos);
                rn1 = reshape(rn1,obj.N,[]);
                rn1 = rn1(obj.tau_order,:);
                rn1 = rn1(obj.keep,:);
                
                rn2 = obj.sector.SpatialConsistency.randomCoupling2{obj.spatialConsist_indx}.rand(obj.UE.initPos);
                rn2 = reshape(rn2,obj.N,[]);
                rn2 = rn2(obj.tau_order,:);
                rn2 = rn2(obj.keep,:);
                
                rn3 = obj.sector.SpatialConsistency.randomCoupling3{obj.spatialConsist_indx}.rand(obj.UE.initPos);
                rn3 = reshape(rn3,obj.N,[]);
                rn3 = rn3(obj.tau_order,:);
                rn3 = rn3(obj.keep,:);
            else
                rn1 = rand(obj.N_new,obj.M);
                rn2 = rand(obj.N_new,obj.M);
                rn3 = rand(obj.N_new,obj.M);
            end
            for n = 1:obj.N_new
                if ~isempty(find(obj.strong_cluster_id == n, 1))
                    % for sub-cluster 1
                    ray_idlist = [1, 2, 3, 4, 5, 6, 7, 8, 19, 20];
%                     [~,index1] = sort(rand(1,10)); 
                    [~,index1] = sort(rn1(n,ray_idlist));
%                         obj.phi_n_m_AOD(n, ray_idlist) = obj.phi_n_m_AOD(n, ray_idlist(index1));
                    obj.phi_n_m_AOA(n, ray_idlist) = obj.phi_n_m_AOA(n, ray_idlist(index1));
%                     [~,index1] = sort(rand(1,10));  
                    [~,index1] = sort(rn2(n,ray_idlist));
%                         obj.theta_n_m_ZOD(n, ray_idlist) = obj.theta_n_m_ZOD(n, ray_idlist(index1));
                    obj.theta_n_m_ZOA(n, ray_idlist) = obj.theta_n_m_ZOA(n, ray_idlist(index1));
%                     [~,index1] = sort(rand(1,10)); 
                    [~,index1] = sort(rn3(n,ray_idlist));
                    obj.phi_n_m_AOD(n, ray_idlist) = obj.phi_n_m_AOD(n, ray_idlist(index1));
%                         obj.theta_n_m_ZOD(n, ray_idlist) = obj.theta_n_m_ZOD(n, ray_idlist(index1));
                    % for sub-cluster 2
                    ray_idlist = [9, 10, 11, 12, 17, 18];
%                     [~,index2] = sort(rand(1,6)); 
                    [~,index2] = sort(rn1(n,ray_idlist));
%                         obj.phi_n_m_AOD(n, ray_idlist) = obj.phi_n_m_AOD(n, ray_idlist(index2));
                    obj.phi_n_m_AOA(n, ray_idlist) = obj.phi_n_m_AOA(n, ray_idlist(index2));
%                     [~,index2] = sort(rand(1,6));  
                    [~,index2] = sort(rn2(n,ray_idlist));
%                         obj.theta_n_m_ZOD(n, ray_idlist) = obj.theta_n_m_ZOD(n, ray_idlist(index2));
                    obj.theta_n_m_ZOA(n, ray_idlist) = obj.theta_n_m_ZOA(n, ray_idlist(index2));
%                     [~,index2] = sort(rand(1,6));  
                    [~,index2] = sort(rn3(n,ray_idlist));
                    obj.phi_n_m_AOD(n, ray_idlist) = obj.phi_n_m_AOD(n, ray_idlist(index2));
%                         obj.theta_n_m_ZOD(n, ray_idlist) = obj.theta_n_m_ZOD(n, ray_idlist(index2));
                    % for sub-cluster 3
                    ray_idlist = [13, 14, 15, 16];
%                     [~,index3] = sort(rand(1,4)); 
                    [~,index3] = sort(rn1(n,ray_idlist));
%                         obj.phi_n_m_AOD(n, ray_idlist) = obj.phi_n_m_AOD(n, ray_idlist(index3));
                    obj.phi_n_m_AOA(n, ray_idlist) = obj.phi_n_m_AOA(n, ray_idlist(index3));
%                     [~,index3] = sort(rand(1,4)); 
                    [~,index3] = sort(rn2(n,ray_idlist));
%                         obj.theta_n_m_ZOD(n, ray_idlist) = obj.theta_n_m_ZOD(n, ray_idlist(index3));
                    obj.theta_n_m_ZOA(n, ray_idlist) = obj.theta_n_m_ZOA(n, ray_idlist(index3));
%                     [~,index3] = sort(rand(1,4)); 
                    [~,index3] = sort(rn3(n,ray_idlist));
                    obj.phi_n_m_AOD(n, ray_idlist) = obj.phi_n_m_AOD(n, ray_idlist(index3));
%                         obj.theta_n_m_ZOD(n, ray_idlist) = obj.theta_n_m_ZOD(n, ray_idlist(index3));
                else
%                     [~,index] = sort(rand(1,obj.M)); 
                    [~,index] = sort(rn1(n,:));
%                         obj.phi_n_m_AOD(n, :) = obj.phi_n_m_AOD(n, index);
                    obj.phi_n_m_AOA(n, :) = obj.phi_n_m_AOA(n, index);
%                     [~,index] = sort(rand(1,obj.M)); 
                    [~,index] = sort(rn2(n,:));
%                         obj.theta_n_m_ZOD(n, :) = obj.theta_n_m_ZOD(n, index);
                    obj.theta_n_m_ZOA(n, :) = obj.theta_n_m_ZOA(n, index);
%                     [~,index] = sort(rand(1,obj.M));
                    [~,index] = sort(rn3(n,:));
                    obj.phi_n_m_AOD(n, :) = obj.phi_n_m_AOD(n, index);
%                         obj.theta_n_m_ZOD(n, :) = obj.theta_n_m_ZOD(n, index);                        
                end
            end
        end
        
       %% calculate the RSRP
%        function [couplingloss,RSRP] = RSRP_cal(obj)
%             tx_power       = obj.sector.Tx_power;
%             M_             = obj.M;
% %            N_             = obj.N_new;
%             UE_            = obj.UE;
%             k_ue           = UE_.antenna.panel(1).K;
%             U              = UE_.antenna.panel(1).num_element/k_ue;
% %             pos_UE         = [UE_.pos, UE_.h_UT]';
%             BS_            = obj.sector;
%             lambda         = 3e8/obj.fc;
% %             pos_BS         = [obj.BS_pos_wrap, BS_.attached_BS.h_BS]';
%             K_             = BS_.antenna.panel(1).K;
%             p              = BS_.antenna.panel(1).num_element/K_;
%             theta_n_m_ZOA_ = obj.theta_n_m_ZOA;
%             theta_n_m_ZOD_ = obj.theta_n_m_ZOD;
%             phi_n_m_AOA_   = obj.phi_n_m_AOA;
%             phi_n_m_AOD_   = obj.phi_n_m_AOD;
%             theta_LOS_ZOA_ = obj.theta_LOS_ZOA;
%             theta_LOS_ZOD_ = obj.theta_LOS_ZOD;
%             phi_LOS_AOA_   = obj.phi_LOS_AOA;
%             phi_LOS_AOD_   = obj.phi_LOS_AOD;
%             PHI_n_m_       = obj.PHI_n_m;
%             PHI_LOS_       = obj.PHI_LOS;
%             XPR_n_m_       = obj.XPR_n_m;
%             H              = zeros([U, p]);  
% 
%             if obj.bLOS && ~obj.O2I
%                 KR = 10.^(obj.K/10);
%             else
%                 KR = 0;
%             end
% 
%             w_ue = UE_.antenna.w_m;
%             r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA_).*cosd(phi_n_m_AOA_); % the spherical unit vector
%             r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA_).*sind(phi_n_m_AOA_);
%             r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA_);
%             
%             
%             w_m  = BS_.antennaArray.w_m;
%             r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD_).*cosd(phi_n_m_AOD_); 
%             r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD_).*sind(phi_n_m_AOD_);
%             r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD_);
%             
%             XPR_nm(1,1,:,:) = squeeze(XPR_n_m_(:,:,1));
% %             XPR_nm2(1,1,:,:) = squeeze(XPR_n_m_(:,:,2));
%             MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m_(1,1,:,:));  MAT_PHI(1,2,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m_(1, 2,:,:));
%             MAT_PHI(2,1,:,:) = sqrt(1./(XPR_nm)).*exp(1j*PHI_n_m_(2, 1,:,:));      MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m_(2, 2,:,:));
%             Pwt = sqrt(obj.Pn/M_/(KR+1));
%             Pwsqrt  = repmat(Pwt',1,M_); 
%                     
%             for u = 1:U  
%                 
%                 kue = (1+k_ue*(u-1)):k_ue*u;
% %                 id = reshape([UE_.antenna.element_list(kue).ID],2,[])';
% %                 d_rx_s = [zeros(length(kue),1), (id(:,2)-1)*UE_.antenna.d_H, (id(:,1)-1)*UE_.antenna.d_V]'*lambda + pos_UE;
%                 d_rx_s = [UE_.antenna.element_list(kue).pos]*lambda;
%                 F_rx_theta_phi_tmp = UE_.antenna.element_list(kue(1)).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi_tmp =  multiprod(multiprod(exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2]))),w_ue.',[1 2],[1 2]),F_rx_theta_phi_tmp,[1 2],[1 2]);
%                 F_rx_theta_phi = permute(F_rx_theta_phi_tmp,[2 1 3 4]);
%                 
%                 for s = 1:p
% 
%                     k = (1+K_*(s-1)):K_*s;
% %                     id = reshape([BS_.antennaArray.element_list(k).ID],2,[])';
% %                     d_tx_s = [zeros(length(k),1), (id(:,2)-1)*BS_.antennaArray.d_H, (id(:,1)-1)*BS_.antennaArray.d_V]'*lambda + pos_BS;
%                     d_tx_s = [BS_.antennaArray.element_list(k).pos]*lambda;
%                     F_tx_theta_phi_tmp = BS_.antennaArray.element_list(k(1)).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                     F_tx_theta_phi =  multiprod(multiprod(exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2]))),w_m.',[1 2],[1 2]),F_tx_theta_phi_tmp,[1 2],[1 2]);
% 
%                     h_tmp = (abs(Pwsqrt.*squeeze(multiprod(multiprod(F_rx_theta_phi,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi,[1 2],[1 2])))).^2;
% 
%                     H(u, s) = sum(sum(h_tmp,2),1);
%                 end
%             end
%             % for LOS case
%             h_LOS = zeros(U,p);
%             if obj.bLOS && ~obj.O2I
%                 r_rx_LOS = [sind(theta_LOS_ZOA_)*cosd(phi_LOS_AOA_);  % the spherical unit vector
%                             sind(theta_LOS_ZOA_)*sind(phi_LOS_AOA_); 
%                             cosd(theta_LOS_ZOA_)];
%                 r_tx_LOS = [sind(theta_LOS_ZOD_)*cosd(phi_LOS_AOD_);  % the spherical unit vector
%                             sind(theta_LOS_ZOD_)*sind(phi_LOS_AOD_); 
%                             cosd(theta_LOS_ZOD_)];
%                 for u = 1:U
%                     kue = (1+k_ue*(u-1)):k_ue*u;
% %                     id = reshape([UE_.antenna.element_list(kue).ID],2,[])';
% %                     d_rx_s = [zeros(length(kue),1), (id(:,2)-1)*UE_.antenna.d_H, (id(:,1)-1)*UE_.antenna.d_V]'*lambda + pos_UE;
%                     d_rx_s = [UE_.antenna.element_list(kue).pos]*lambda;
%                     F_rx_theta_phi = UE_.antenna.element_list(kue(1)).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_rx_theta_phi = sum((w_ue.*exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s))),2)*F_rx_theta_phi;
%                     for s =1:p
%                         k = (1+K_*(s-1)):K_*s;
% %                         id = reshape([BS_.antennaArray.element_list(k).ID],2,[])';
% %                         d_tx_s = [zeros(length(k),1), (id(:,2)-1)*BS_.antennaArray.d_H, (id(:,1)-1)*BS_.antennaArray.d_V]'*lambda + pos_BS;
%                         d_tx_s = [BS_.antennaArray.element_list(k).pos]*lambda;
%                         F_tx_theta_phi_tmp = BS_.antennaArray.element_list(k(1)).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                         F_tx_theta_phi = sum((w_m.*exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s))),2)*F_tx_theta_phi_tmp;
%                         h_LOS(u,s) = (abs(sqrt(KR/(KR+1))*(F_rx_theta_phi.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi))).^2; 
%                     end
%                 end
%             end
%             couplingloss = -obj.PL + obj.SF + 10*log10(sum((h_LOS + H),1)) - 10*log10(U);% + tx_power;
%             RSRP         = -obj.PL + obj.SF + 10*log10(sum((h_LOS + H),1)) - 10*log10(U) + tx_power;
%        end
%        
       
       %% Step 11: Generate channel coefficients for each cluster n and each receiver and transmitter element pair u, s.
%        function [H, tau] = channel_coef_gen(obj, t, largeScaleEnable)
%             M_             = obj.M;
%             N_             = obj.N_new;
%             UE_            = obj.UE;
% %             pos_UE         = [UE_.pos, UE_.h_UT]';
%             U              = UE_.antenna.element_num;
%             BS_            = obj.sector;
% %             pos_BS         = [obj.BS_pos_wrap, BS_.attached_BS.h_BS]';
%             S              = BS_.antennaArray.element_num;
%             theta_n_m_ZOA_ = obj.theta_n_m_ZOA;
%             theta_n_m_ZOD_ = obj.theta_n_m_ZOD;
%             phi_n_m_AOA_   = obj.phi_n_m_AOA;
%             phi_n_m_AOD_   = obj.phi_n_m_AOD;
%             theta_LOS_ZOA_ = obj.theta_LOS_ZOA;
%             theta_LOS_ZOD_ = obj.theta_LOS_ZOD;
%             phi_LOS_AOA_   = obj.phi_LOS_AOA;
%             phi_LOS_AOD_   = obj.phi_LOS_AOD;
%             PHI_n_m_       = obj.PHI_n_m;
%             PHI_LOS_       = obj.PHI_LOS;
%             XPR_n_m_       = obj.XPR_n_m;
% 
%             lambda         = 3e8/obj.fc;
%             % N+4: (the N¨C2 weakest clusters and the two strongest clusters with three sub-cluster) 
%             tau            = zeros(1,N_+4);
% 
%             v_vec          = UE_.v*[sind(UE_.theta_v)*cosd(UE_.phi_v);
%                                   sind(UE_.theta_v)*sind(UE_.phi_v);
%                                   cosd(UE_.theta_v)];
% 
%             if obj.bLOS && ~obj.O2I
%                 tau_tmp = obj.tau_n_LOS;
%                 KR = 10.^(obj.K/10);
%             else
%                 tau_tmp = obj.tau_n;
%                 KR = 0;
%             end
%             XPR_nm(1,1,:,:)  = squeeze(XPR_n_m_(:,:,1));
% %             XPR_nm2(1,1,:,:) = squeeze(XPR_n_m_(:,:,2));
%             MAT_PHI(1,1,:,:) = exp(1j*PHI_n_m_(1,1,:,:));                  MAT_PHI(1,2,:,:) = sqrt(1./XPR_nm).*exp(1j*PHI_n_m_(1, 2,:,:));
%             MAT_PHI(2,1,:,:) = sqrt(1./XPR_nm).*exp(1j*PHI_n_m_(2, 1,:,:)); MAT_PHI(2,2,:,:) = exp(1j*PHI_n_m_(2, 2,:,:));
%             Pwt              = sqrt(obj.Pn/M_/(KR+1));
%             Pwsqrt           = repmat(Pwt',1,M_);
%             
%             map_choise = zeros([U,S,N_,M_,N_+4]);
%             path_id       = 0;
%             for n = 1:N_
%                 if ~isempty(find(obj.strong_cluster_id == n, 1))
%                     iter_num = 3;
%                 else
%                     iter_num = 1;
%                 end
%                 for i = 1:iter_num
%                     % for next path
%                     path_id = path_id + 1;
%                     
%                     if i == 1
%                         tau(path_id) = tau_tmp(n);
%                         if ~isempty(find(obj.strong_cluster_id == n, 1))
%                             m_list = [1, 2, 3, 4, 5, 6, 7, 8, 19, 20];
%                         else
%                             m_list = 1:M_;
%                         end
%                     elseif i == 2
%                         m_list = [9, 10, 11, 12, 17, 18];
%                         tau(path_id) = tau_tmp(n)+5e-9;
%                     else
%                         m_list = [13, 14, 15, 16];
%                         tau(path_id) = tau_tmp(n)+10e-9;
%                     end
%                     map_choise(:,:,n,m_list,path_id) = 1;
%                 end
%             end
%             H = zeros([U, S, path_id]);
%             tau = tau(1:path_id);
% 
% %                 w_ue = UE_.ant.w_m;
%             r_rx_n_m(1,1,:,:) = sind(theta_n_m_ZOA_).*cosd(phi_n_m_AOA_); % the spherical unit vector
%             r_rx_n_m(1,2,:,:) = sind(theta_n_m_ZOA_).*sind(phi_n_m_AOA_);
%             r_rx_n_m(1,3,:,:) = cosd(theta_n_m_ZOA_);
% %                 w_m  = BS_.ant.w_m;
%             r_tx_n_m(1,1,:,:) = sind(theta_n_m_ZOD_).*cosd(phi_n_m_AOD_); 
%             r_tx_n_m(1,2,:,:) = sind(theta_n_m_ZOD_).*sind(phi_n_m_AOD_);
%             r_tx_n_m(1,3,:,:) = cosd(theta_n_m_ZOD_);
%             
% %             id = reshape([UE_.antenna.element_list(1:U).ID],2,[])';
% %             d_rx_s = [zeros([U,1]), (id(:,2)-1)*UE_.antenna.d_H, (id(:,1)-1)*UE_.antenna.d_V]'*lambda + pos_UE;
%             d_rx_s = [UE_.antenna.element_list(1:U).pos]*lambda;
% %             id = reshape([BS_.antennaArray.element_list(1:S).ID],2,[])';
% %             d_tx_s = [zeros([S,1]), (id(:,2)-1)*BS_.antennaArray.d_H, (id(:,1)-1)*BS_.antennaArray.d_V]'*lambda + pos_BS;
%             d_tx_s = [BS_.antennaArray.element_list(1:S).pos]*lambda;
%             if UE_.antenna.P == 1 && BS_.antennaArray.P == 1
%                 F_rx_theta_phi_tmp = UE_.antenna.element_list(1).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi1    = permute(F_rx_theta_phi_tmp,[2 1 3 4]);
%                 F_tx_theta_phi1    = BS_.antennaArray.element_list(1).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                 F_total            = squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi1,[1 2],[1 2]));
%                 F_total_e          = permute(repmat(F_total,1,1,U,S),[3,4,1,2]);
%             elseif UE_.antenna.P == 2 && BS_.antennaArray.P == 1
%                 F_rx_theta_phi_tmp1 = UE_.antenna.element_list(1).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
%                 F_rx_theta_phi_tmp2 = UE_.antenna.element_list(U).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
%                 F_tx_theta_phi1     = BS_.antennaArray.element_list(1).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                 F_total_1           = squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi1,[1 2],[1 2]));
%                 F_total_2           = squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi1,[1 2],[1 2]));
%                 F_total_e           = permute(cat(3,repmat(F_total_1,1,1,U/2,S),repmat(F_total_2,1,1,U/2,S)),[3,4,1,2]);
%                 
%             elseif UE_.antenna.P == 1 && BS_.antennaArray.P == 2
%                 F_rx_theta_phi_tmp1 = UE_.antenna.element_list(1).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
%                 F_tx_theta_phi1     = BS_.antennaArray.element_list(1).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                 F_tx_theta_phi2     = BS_.antennaArray.element_list(S).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                 F_total_1           = squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi1,[1 2],[1 2]));
%                 F_total_2           = squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi2,[1 2],[1 2]));
%                 F_total_e           = permute(cat(4,repmat(F_total_1,1,1,U,S/2),repmat(F_total_2,1,1,U,S/2)),[3,4,1,2]);
%                 
%             elseif UE_.antenna.P == 2 && BS_.antennaArray.P == 2
%                 F_rx_theta_phi_tmp1 = UE_.antenna.element_list(1).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi1     = permute(F_rx_theta_phi_tmp1,[2 1 3 4]);
%                 F_rx_theta_phi_tmp2 = UE_.antenna.element_list(U).field_pattern(phi_n_m_AOA_, theta_n_m_ZOA_);
%                 F_rx_theta_phi2     = permute(F_rx_theta_phi_tmp2,[2 1 3 4]);
%                 F_tx_theta_phi1     = BS_.antennaArray.element_list(1).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                 F_tx_theta_phi2     = BS_.antennaArray.element_list(S).field_pattern(phi_n_m_AOD_, theta_n_m_ZOD_);
%                 F_total_1           = squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi1,[1 2],[1 2]));
%                 F_total_2           = squeeze(multiprod(multiprod(F_rx_theta_phi1,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi2,[1 2],[1 2]));
%                 F_total_3           = squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi1,[1 2],[1 2]));
%                 F_total_4           = squeeze(multiprod(multiprod(F_rx_theta_phi2,MAT_PHI,[1 2],[1 2]),F_tx_theta_phi2,[1 2],[1 2]));
%                 F_total_12          = cat(4,repmat(F_total_1,1,1,U/2,S/2),repmat(F_total_2,1,1,U/2,S/2));
%                 F_total_34          = cat(4,repmat(F_total_3,1,1,U/2,S/2),repmat(F_total_4,1,1,U/2,S/2));
%                 F_total_e           = permute(cat(3,F_total_12,F_total_34),[3,4,1,2]);
%                 
%             end
% 
%             phase_rx    = exp(1j*2*pi/lambda*(multiprod(r_rx_n_m,d_rx_s,[1 2],[1 2])));
%             phase_rx_e  = repmat(permute(phase_rx,[2 1 3 4]),1,S,1,1);
%             phase_tx    = exp(1j*2*pi/lambda*(multiprod(r_tx_n_m,d_tx_s,[1 2],[1 2])));
%             phase_tx_e  = repmat(phase_tx,U,1,1,1);
%             phase_vnm   = exp(1j*2*pi/lambda*multiprod(r_rx_n_m,v_vec,[1,2],[1,2])*t);
%             phase_vnm_e = repmat(phase_vnm,U,S,1,1);
% 
%             F_e = F_total_e.*phase_rx_e.*phase_tx_e.*phase_vnm_e;
%             H_e = zeros([U, S, N_, M_]);
% 
%             for u = 1:U
%                 for s = 1:S
%                     H_e(u,s,:,:) = squeeze(F_e(u,s,:,:)).*Pwsqrt;
%                 end
%             end
%             for pa = 1:path_id
%                 H(:,:,pa) = squeeze(sum(sum(H_e.*squeeze(map_choise(:,:,:,:,pa)),4),3));
%             end
% 
%             % for LOS case
%             if obj.bLOS && ~obj.O2I
%                 r_rx_LOS = [sind(theta_LOS_ZOA_)*cosd(phi_LOS_AOA_);  % the spherical unit vector
%                             sind(theta_LOS_ZOA_)*sind(phi_LOS_AOA_); 
%                             cosd(theta_LOS_ZOA_)];
%                 r_tx_LOS = [sind(theta_LOS_ZOD_)*cosd(phi_LOS_AOD_);  % the spherical unit vector
%                             sind(theta_LOS_ZOD_)*sind(phi_LOS_AOD_); 
%                             cosd(theta_LOS_ZOD_)];
% 
% %                 id = reshape([UE_.antenna.element_list(1:U).id],2,[])';
% %                 d_rx_s = [zeros([U,1]), (id(:,2)-1)*UE_.antenna.d_H, (id(:,1)-1)*UE_.antenna.d_V]'*lambda + pos_UE;
%                 d_rx_s = [UE_.antenna.element_list(1:U).pos]*lambda;
% %                 id = reshape([BS_.antennaArray.element_list(1:S).id],2,[])';
% %                 d_tx_s = [zeros([S,1]), (id(:,2)-1)*BS_.antennaArray.d_H, (id(:,1)-1)*BS_.antennaArray.d_V]'*lambda + pos_BS;
%                 d_tx_s = [BS_.antennaArray.element_list(1:S).pos]*lambda;
%                 if UE_.antenna.P == 1 && BS_.antennaArray.P == 1
%                     
%                     F_rx_theta_phi1 = UE_.antenna.element_list(1).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_tx_theta_phi1 = BS_.antennaArray.element_list(1).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                     F_total         = sqrt(KR/(KR+1))*F_rx_theta_phi1.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi1;
%                     F_total_e       = repmat(F_total,U,S);
%                 elseif UE_.antenna.P == 2 && BS_.antennaArray.P == 1
%                     
%                     F_rx_theta_phi1 = UE_.antenna.element_list(1).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_rx_theta_phi2 = UE_.antenna.element_list(U).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_tx_theta_phi1 = BS_.antennaArray.element_list(1).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                     F_total_1       = sqrt(KR/(KR+1))*F_rx_theta_phi1.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi1;
%                     F_total_2       = sqrt(KR/(KR+1))*F_rx_theta_phi2.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi1;
%                     F_total_e       = [repmat(F_total_1,U/2,S);repmat(F_total_2,U/2,S)];
%                 elseif UE_.antenna.P == 1 && BS_.antennaArray.P == 2
%                     
%                     F_rx_theta_phi1 = UE_.antenna.element_list(1).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_tx_theta_phi1 = BS_.antennaArray.element_list(1).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                     F_tx_theta_phi2 = BS_.antennaArray.element_list(S).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                     F_total_1       = sqrt(KR/(KR+1))*F_rx_theta_phi1.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi1;
%                     F_total_2       = sqrt(KR/(KR+1))*F_rx_theta_phi1.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi2;
%                     F_total_e       = [repmat(F_total_1,U,S/2),repmat(F_total_2,U,S/2)];
%                 elseif UE_.antenna.P == 2 && BS_.antennaArray.P == 2
%                     
%                     F_rx_theta_phi1 = UE_.antenna.element_list(1).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_rx_theta_phi2 = UE_.antenna.element_list(U).field_pattern(phi_LOS_AOA_, theta_LOS_ZOA_);
%                     F_tx_theta_phi1 = BS_.antennaArray.element_list(1).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                     F_tx_theta_phi2 = BS_.antennaArray.element_list(S).field_pattern(phi_LOS_AOD_, theta_LOS_ZOD_);
%                     F_total_1       = sqrt(KR/(KR+1))*F_rx_theta_phi1.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi1;
%                     F_total_2       = sqrt(KR/(KR+1))*F_rx_theta_phi1.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi2;
%                     F_total_3       = sqrt(KR/(KR+1))*F_rx_theta_phi2.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi1;
%                     F_total_4       = sqrt(KR/(KR+1))*F_rx_theta_phi2.'*[exp(1j*PHI_LOS_(1)), 0; 0, -exp(1j*PHI_LOS_(1))]*F_tx_theta_phi2;
%                     F_total_e       = [[repmat(F_total_1,U/2,S/2),repmat(F_total_2,U/2,S/2)];[repmat(F_total_3,U/2,S/2),repmat(F_total_4,U/2,S/2)]];
%                     
%                 end
% 
%                 phase_rx    = exp(1j*2*pi/lambda*(r_rx_LOS'*d_rx_s));
%                 phase_rx_e  = repmat(phase_rx.',1,S);
%                 phase_tx    = exp(1j*2*pi/lambda*(r_tx_LOS'*d_tx_s));
%                 phase_tx_e  = repmat(phase_tx,U,1);
%                 phase_vnm   = exp(1j*2*pi/lambda*(r_rx_LOS'*v_vec)*t);
%                 phase_vnm_e = repmat(phase_vnm,U,S);
%                 H(:,:,1)    = H(:,:,1) + F_total_e.*phase_rx_e.*phase_tx_e.*phase_vnm_e;
%             end
%                 
%            %% Step 12: Apply pathloss and shadowing for the channel coefficients
%             if largeScaleEnable
%                 H = H*(sqrt(10.^(-obj.PL/10))*sqrt(10.^(obj.SF/10)));
%             end
%        end
%        
       function [ as, mean_angle ]  = calc_angular_spreads(obj,ang,pow,LOSang,wrap_angles )
            % CALC_ANGULAR_SPREADS Calculates the angular spread in degree

            ang = ang(:).'*pi/180;
            pow = repmat((pow/obj.M).',1,obj.M).*obj.loss_blockage;
            pow = pow(:).';
            % Normalize powers
            pt = sum( pow,2 );
            pow = pow./pt;
            if obj.bLOS && ~obj.O2I
                KF = 10^(obj.K/10);
                pow = pow*(1/(1+KF));
                ang(numel(pow)+1) = LOSang*pi/180;
                pow(numel(pow)+1) = KF/(1+KF);
            end
            if ~exist('wrap_angles','var')
                wrap_angles = true;
            end
            if wrap_angles
                mean_angle = angle( sum( pow.*exp( 1j*ang ) , 2 ) ); % [rad]
            else
                mean_angle = sum( pow.*ang,2 ); 
            end
            phi = ang - mean_angle;
            if wrap_angles
                phi = angle( exp( 1j*phi ) );
            end
%             as1 = sqrt(-2*log(abs(sum( pow.*exp( 1j*ang ),2))));
            as = sqrt( sum(pow.*(phi.^2),2)); %  - sum( pow.*phi,2).^2 
%             as  = min(as1,as2);
            mean_angle = mean_angle*180/pi;
            as = as*180/pi;
       end
       
       function [ ds, mean_delay ] = calc_delay_spread(obj, taus, pow )
       % CALC_DELAY_SPREAD Calculates the delay spread in [s]
            
            taus = repmat( taus',1,obj.M ); %  + obj.d_3D/3e8
            pow = (obj.loss_blockage(:,1).*pow(:)).';
            % Normalize powers
            pt = sum( pow,2 );
            pow = pow./pt;
%             if obj.bLOS && ~obj.O2I
%                 KF = 10^(obj.K/10);
%                 pow = pow*(1/(1+KF));
%                 pow(1) = KF/(1+KF)+pow(1);
%             end
            pow = repmat( (pow/obj.M).',1,obj.M );
            taus(obj.strong_cluster_id,:) = taus(obj.strong_cluster_id,:)+repmat(obj.map_delay,numel(obj.strong_cluster_id),1);
            pow = pow(:);
            taus = taus(:);
            if ~isempty(obj.oxygenAbsorp)
                [~,loss_oxy] = obj.oxygenAbsorp.loss(obj.d_3D, taus);
            else
                loss_oxy = 1;
            end
            pow = loss_oxy.*pow;
            pow = pow/sum(pow);
            
            mean_delay = sum( pow.*taus ); 
            tmp = taus - mean_delay;
%             ds = sqrt( sum(pow.*((tmp).^2)) ); %  - sum(pow.*tmp).^2
            ds = sqrt( sum(pow.*((taus).^2)) - mean_delay.^2 ); %  - sum(pow.*tmp).^2 
        end
    end
    
end

% other function
function [phi, theta] = cart2sph_(x, y, z)
    [phi,ele,~] = cart2sph(x,y,z);
    phi         = phi*180/pi;          % [-180, 180]
    theta       = 90 - ele*180/pi;     % [   0, 180]
end
function [xyz] = sph2cart_(phi, theta, varargin)
    if isempty(varargin)
        r = ones(size(phi));
    else 
        r = varargin{1};
    end
    x              = r.*sind(theta).*cosd(phi);
    y              = r.*sind(theta).*sind(phi);
    z              = r.*cosd(theta);
    [n, m]         = size(phi);
    xyz(1,1:n,1:m) = x;
    xyz(2,1:n,1:m) = y;
    xyz(3,1:n,1:m) = z;
%     xyz = [x y z];
end