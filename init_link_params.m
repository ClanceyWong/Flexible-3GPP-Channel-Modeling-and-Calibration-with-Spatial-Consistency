function link_params = init_link_params(obj)
    fc = obj.frequency;
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
                obj.mu_offset_ZOD = (7.66*log10(fc)-5.96) -10^((0.208*log10(fc)-0.782)*log10(max((25), obj.d_2D_out))+(-0.13*log10(fc)+2.03)-0.07*(obj.UE.h_UT-1.5));
            end
            temp1 = obj.sector.attached_BS.corr_sos{3}.randn([obj.UE.pos,obj.UE.h_UT]'); % obj.UE.n_fl
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.UMa_C_sqrt_O2I*ksi;
            % SF/DS/ASD/ASA/ZSD/ZSA
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
            obj.c_ASA = 8;
            obj.c_ZSA = 3;
            obj.c_ASD = 5;
        elseif obj.bLOS % 0 1
            obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.UE.h_UT-1.5)+0.75);
            obj.sigma_lgZSD = 0.4;
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
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
            obj.mu_offset_ZOD = (7.66*log10(fc)-5.96) -10^((0.208*log10(fc)-0.782)*log10(max((25), obj.d_2D))+(-0.13*log10(fc)+2.03)-0.07*(obj.UE.h_UT-1.5));
            temp1 = obj.sector.attached_BS.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
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
                obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D_out/1000)+0.01*abs(obj.UE.h_UT-obj.sector.attached_BS.h_BS)+0.83);
                obj.sigma_lgZSD = 0.35; % 0.4
                obj.mu_offset_ZOD = 0;
            else
                % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*max(obj.UE.h_UT-obj.sector.attached_BS.h_BS, 0)+0.9);
                obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D_out/1000)+0.01*max(obj.UE.h_UT-obj.sector.attached_BS.h_BS, 0)+0.2);
                obj.sigma_lgZSD = 0.35; %0.6
                % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D_out))+1.6);
                obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D_out))+3.3);
            end
            temp1 = obj.sector.attached_BS.corr_sos{3}.randn([obj.UE.pos,obj.UE.h_UT]');
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
            obj.c_ASA = 8;
            obj.c_ZSA = 3;
            obj.c_ASD = 5;
        elseif obj.bLOS % 0 1
            % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*abs(obj.UE.h_UT-obj.sector.attached_BS.h_BS)+0.75);
            obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D/1000)+0.01*abs(obj.UE.h_UT-obj.sector.attached_BS.h_BS)+0.83);
            obj.sigma_lgZSD = 0.35; % 0.4
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
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
            obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D/1000)+0.01*max(obj.UE.h_UT-obj.sector.attached_BS.h_BS, 0)+0.2);
            obj.sigma_lgZSD = 0.6;
            % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D))+1.6);
            obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D))+3.3);
            temp1 = obj.sector.attached_BS.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
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
    elseif strcmp(obj.scenario.name,'3D-RMa')
        if obj.O2I
            obj.mu_lgZSD = max(-1, -0.19*(obj.d_2D_out/1000)-0.01*(obj.UE.h_UT-1.5)+0.28);
            obj.sigma_lgZSD = 0.30;
            obj.mu_offset_ZOD = atand((35-3.5)/obj.d_2D_out)-atand((35-1.5)/obj.d_2D_out);
            temp1 = obj.sector.attached_BS.corr_sos{3}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.RMa_C_sqrt_O2I*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.DS = 10^(0.24*ksi(2)-7.47);
            obj.ASD = 10^(0.18*ksi(3)+0.67);
            obj.ASA = 10^(0.21*ksi(4)+1.66);
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
            obj.ZSA = 10^(0.22*ksi(6)+0.93);
            obj.r_tau = 1.7;
            obj.mu_XPR = 7;
            obj.sigma_XPR = 3;
            obj.N = 10;
            obj.M = 20;
            obj.zeta = 3;
            obj.c_ASA = 3;
            obj.c_ZSA = 3;
            obj.c_ASD = 2;
        elseif obj.bLOS
            obj.mu_lgZSD = max(-1, -0.17*(obj.d_2D/1000)-0.01*(obj.UE.h_UT-1.5)+0.22);
            obj.sigma_lgZSD = 0.34;
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.RMa_C_sqrt_LOS*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.K = 4*ksi(2)+7;
            obj.DS = 10^(0.55*ksi(3)-7.49);
            obj.ASD = 10^(0.38*ksi(4)+0.90);
            obj.ASA = 10^(0.24*ksi(5)+1.52);
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
            obj.ZSA = 10^(0.40*ksi(7)+0.47);
            obj.r_tau = 3.8;
            obj.mu_XPR = 12;
            obj.sigma_XPR = 4;
            obj.N = 11;
            obj.M = 20;
            obj.zeta = 3;
            obj.c_ASA = 3;
            obj.c_ZSA = 2;
            obj.c_ASD = 2;
        else
            obj.mu_lgZSD = max(-1, -0.19*(obj.d_2D/1000)-0.01*(obj.UE.h_UT-1.5)+0.28);
            obj.sigma_lgZSD = 0.30;
            obj.mu_offset_ZOD = atand((35-3.5)/obj.d_2D)-atand((35-1.5)/obj.d_2D);
            temp1 = obj.sector.attached_BS.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.RMa_C_sqrt_NLOS*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.DS = 10^(0.48*ksi(2)-7.43);
            obj.ASD = 10^(0.45*ksi(3)+0.95);
            obj.ASA = 10^(0.13*ksi(4)+1.52);
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
            obj.ZSA = 10^(0.37*ksi(6)+0.58);
            obj.r_tau = 1.7;
            obj.mu_XPR = 7;
            obj.sigma_XPR = 3;
            obj.N = 10;
            obj.M = 20;
            obj.zeta = 3;
            obj.c_ASA = 3;
            obj.c_ZSA = 2.5;
            obj.c_ASD = 2;
        end
    elseif  strcmp(obj.scenario.name,'3D-InH')
        fc(fc<=6) = 6;
        if obj.bLOS
            % obj.mu_lgZSD = 1.02;
            % obj.sigma_lgZSD = 0.41;
            obj.mu_lgZSD = -1.43*log10(1+fc)+2.228;
            obj.sigma_lgZSD = 0.13*log10(1+fc)+0.3;
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.InH_C_sqrt_LOS*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.K = 4*ksi(2)+7;
            obj.DS = 10^(0.18*ksi(3)+ (-0.01*log10(1+fc)-7.692));
            obj.ASD = 10^(0.18*ksi(4)+1.60);
            obj.ASA = 10^((0.12*log10(1+fc)+0.119)*ksi(5)+ (-0.19*log10(1+fc)+1.781));
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
            obj.ZSA = 10^((-0.04*log10(1+fc)+0.264)*ksi(7)+ (-0.26*log10(1+fc)+1.44));
            obj.r_tau = 3.6;
            obj.mu_XPR = 11;
            obj.sigma_XPR = 4; % 3
            obj.N = 15;
            obj.M = 20;
            obj.zeta = 6;
            obj.c_ASA = 8;
            obj.c_ZSA = 9;
            obj.c_ASD = 5;
        else
            obj.mu_lgZSD = 1.08;
            obj.sigma_lgZSD = 0.36;
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.InH_C_sqrt_NLOS*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.DS = 10^((0.1*log10(1+fc)+0.055)*ksi(2)+ (-0.28*log10(1+fc)-7.173));
            obj.ASD = 10^(0.25*ksi(3)+1.62);
            obj.ASA = 10^((0.12*log10(1+fc)+0.059)*ksi(4)+ (-0.11*log10(1+fc)+1.863));
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
            obj.ZSA = 10^((-0.09*log10(1+fc)+0.746)*ksi(6)+ (-0.15*log10(1+fc)+1.387));
            obj.r_tau = 3;
            obj.mu_XPR = 10;
            obj.sigma_XPR = 4; % 3
            obj.N = 19;
            obj.M = 20;
            obj.zeta = 3;
            obj.c_ASA = 11;
            obj.c_ZSA = 9;
            obj.c_ASD = 5;
        end
    elseif  strcmp(obj.scenario.name,'InF')
        % V = hall volume in m^3, S = total surface area of hall in m^2 (walls+floor+ceiling)
        if obj.bLOS
            obj.mu_lgZSD = 1.35;
            obj.sigma_lgZSD = 0.35;
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{1}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.InH_C_sqrt_LOS*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.K = 8*ksi(2)+7;
            obj.DS = 10^(0.15*ksi(3)+ (log10(26*(obj.scenario.V/obj.scenario.S)+14)-9.35));
            obj.ASD = 10^(0.25*ksi(4)+1.56);
            obj.ASA = 10^((0.12*log10(1+fc)+0.2)*ksi(5)+ (-0.18*log10(1+fc)+1.78));
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(6)+obj.mu_lgZSD);
            obj.ZSA = 10^(0.35*ksi(7)+ (-0.2*log10(1+fc)+1.5));
            obj.r_tau = 2.7;
            obj.mu_XPR = 12;
            obj.sigma_XPR = 6;
            obj.N = 25;
            obj.M = 20;
            obj.zeta = 4;
            obj.c_ASA = 8;
            obj.c_ZSA = 9;
            obj.c_ASD = 5;
        else
            obj.mu_lgZSD = 1.2;
            obj.sigma_lgZSD = 0.55;
            obj.mu_offset_ZOD = 0;
            temp1 = obj.sector.attached_BS.corr_sos{2}.randn([obj.UE.pos,obj.UE.h_UT]');
            ksi = reshape(temp1,[numel(temp1), 1]);
            ksi = tab.InH_C_sqrt_NLOS*ksi;
            obj.SF = obj.sigma_SF*ksi(1);
            obj.DS = 10^(0.19*ksi(2)+ (log10(30*(obj.scenario.V/obj.scenario.S)+32)-9.44));
            obj.ASD = 10^(0.2*ksi(3)+1.57);
            obj.ASA = 10^(0.3*ksi(4)+1.72);
            obj.ZSD = 10^(obj.sigma_lgZSD*ksi(5)+obj.mu_lgZSD);
            obj.ZSA = 10^(0.45*ksi(6)+ (-0.13*log10(1+fc)+1.45));
            obj.r_tau = 3;
            obj.mu_XPR = 11;
            obj.sigma_XPR = 6;
            obj.N = 25;
            obj.M = 20;
            obj.zeta = 3;
            obj.c_ASA = 8;
            obj.c_ZSA = 9;
            obj.c_ASD = 5;
        end
    end

end