plot_flag   = 1;
set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

N0_dBm      = -174;
N0_linear   = (10^(N0_dBm/10));  % mW/Hz 

layer_num   = 3;              % the number of layer of cells.
sec_per_BS  = 3;              % the number of sector per BS.
ue_per_sec  = 30;             % the number of ue per sector.
% load('.\+data\parameters_tab.mat', 'tab');
tab = data.parameters_tab();
UEant_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
UEant_params.array.Ng           = 1;
UEant_params.array.dg_H         = 2.5;
UEant_params.array.dg_V         = 2.5;
UEant_params.panel.M            = 1;                 %  the number of antenna elements with the same polarization in each column.
UEant_params.panel.N            = 1;
UEant_params.panel.Kv           = 1;                  %  K = M or 1
UEant_params.panel.Kh           = 1; 
UEant_params.panel.d_H          = 0.5; 
UEant_params.panel.d_V          = 0.5; 
UEant_params.panel.P            = 2;
UEant_params.panel.X_pol        = [0, 90]; 
UEant_params.panel.ele_downtilt = 90;
UEant_params.panel.ele_panning  = 0;
UEant_params.panel.pol_model    = 'model-1';
beta                            = 0;


today = char(datetime('today'));
% today = '2020-07-17';
if ~exist(['.\result\Phase2_36873_result','_',today],'file') 
    mkdir(['.\result\Phase2_36873_result','_',today]);
end 

bandwidth   = 1e7;           % 10 MHz
for frequency = 2e9
    for case_num = 1:4
        if case_num == 1 || case_num == 2
            BSant_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
            BSant_params.array.Ng           = 1;
            BSant_params.array.dg_H         = 2.5;
            BSant_params.array.dg_V         = 2.5;
            BSant_params.panel.M            = 2;                 %  the number of antenna elements with the same polarization in each column.
            BSant_params.panel.Kv           = 1;                  %  K = M or 1
            BSant_params.panel.Kh           = 1;
            BSant_params.panel.N            = 2;
            BSant_params.panel.d_H          = 0.5; 
            BSant_params.panel.d_V          = 0.5; 
            BSant_params.panel.P            = 1;
            BSant_params.panel.X_pol        = 0; 
            BSant_params.panel.ele_downtilt = 102;
            BSant_params.panel.ele_panning  = 0;
            BSant_params.panel.pol_model    = 'model-1';
            
            UEant_params.panel.N            = 2;
            UEant_params.panel.P            = 1;
            UEant_params.panel.X_pol        = 0; 
        else
            BSant_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
            BSant_params.array.Ng           = 1;
            BSant_params.array.dg_H         = 2.5;
            BSant_params.array.dg_V         = 2.5;
            BSant_params.panel.M            = 10;                 %  the number of antenna elements with the same polarization in each column.
            BSant_params.panel.Kv           = 10;                  %  K = M or 1
            BSant_params.panel.Kh           = 1;
            BSant_params.panel.N            = 2;
            BSant_params.panel.d_H          = 0.5; 
            BSant_params.panel.d_V          = 0.5; 
            BSant_params.panel.P            = 2;
            BSant_params.panel.X_pol        = [-45,45]; 
            BSant_params.panel.ele_downtilt = 102;
            BSant_params.panel.ele_panning  = 0;
            BSant_params.panel.pol_model    = 'model-1';
            
            UEant_params.panel.N            = 1;
            UEant_params.panel.P            = 2;
            UEant_params.panel.X_pol        = [0,90]; 
            beta                            = 90;
        end

        if case_num == 1 || case_num == 3
            scenario = scenarios.UMa_3D(layer_num);
            scenario.BW = bandwidth;
        elseif case_num == 2 || case_num == 4
            scenario = scenarios.UMi_3D(layer_num);
            scenario.BW = bandwidth;
        end

        Noise_mW    = 10^((10*log10(N0_linear * bandwidth) + 9)/10);

        % drop BS
        [BS, BSsector_list]    = network_layout.Drop_BaseStation(scenario, BSant_params, sec_per_BS, frequency, plot_flag);
        % drop UE
        if strcmp(scenario.name,'3D-RMa') ||strcmp(scenario.name,'3D-UMa')||strcmp(scenario.name,'3D-UMi')
            [UE_list, ue_pos_list] = network_layout.Drop_UE(BSsector_list, UEant_params, ue_per_sec, scenario.R, scenario.min_d, beta, plot_flag);
        elseif strcmp(scenario.name,'3D-InH')
            [UE_list, ue_pos_list] = network_layout.Drop_UE_InH(BSsector_list, UEant_params, ue_per_sec, scenario.min_d, [scenario.x_range;scenario.y_range], beta, plot_flag);
        end
        [~, indx]              = sort(ue_pos_list(:,3));

        %% create all links
        port0 = 0;
        t = 0;
        CouplingLoss = zeros(1,numel(UE_list)); 
        RSRP         = zeros(1,numel(UE_list)); 
        Geometry     = zeros(1,numel(UE_list));
        ZSA_ZSD      = zeros(2,numel(UE_list));
        SV           = zeros(2,numel(UE_list));
        RSV          = zeros(1,numel(UE_list));
        wb           = waitbar(0,'creating all links...');
        for nn = 1:numel(UE_list)
            n1        = indx(nn);
%             if UE_list(n1).h_UT ~= UE_list(indx(max(nn-1,1))).h_UT
%                 for fl = 1:numel(BS)
%                     BS(fl).getRandnGrid();
%                 end
%             end
            link_list = [];
            Ue           = UE_list(n1);
            couplingloss = zeros(size(BSsector_list));
            RSRPt        = zeros(size(BSsector_list));
            for n2 = 1:numel(BS)
                site      = BS(n2);
                link      = links.Link(site, Ue, scenario, true, tab);
                for n3 = 1:sec_per_BS
                    [RSRPt(3*(n2-1)+n3), couplingloss(3*(n2-1)+n3)] = Ue.RSRP_calc(BSsector_list(site.sector(n3)),link,port0);
                end
                link_list = [link_list;link]; %#ok<AGROW>
%                 RSRP      = link.RSRP_cal();
%                 PL(n2)    = -link.PL + link.SF + RSRP{1}(1);
%                 link.CouplingLoss = PL(n2);
            end
            [pl, link_id]    = max(RSRPt);
            link             = link_list(ceil(link_id/3));
            
            BSsector_list(link_id).init_link_params(link);
            BSsector_list(link_id).link_params{length(BSsector_list(link_id).link_params)}.PHI_n_m = BSsector_list(link_id).PHI_n_m;
            BSsector_list(link_id).PHI_n_m = [];
            
            RSRP(n1)         = RSRPt(link_id); 
            CouplingLoss(n1) = couplingloss(link_id); 
            RSRPt(link_id)   = -inf;
            interf           = 10*log10(sum(10.^(RSRPt/10)));
            Geometry(n1)     = RSRP(n1)-interf;
            Ue.couplingloss  = CouplingLoss(n1);
            Ue.SIR_geometry  = Geometry(n1);
            
            ZSA_ZSD(:,n1)    = [link.calc_angular_spreads(link.theta_n_m_ZOA,link.Pn,link.theta_LOS_ZOA); ...
                                link.calc_angular_spreads(link.theta_n_m_ZOD,link.Pn,link.theta_LOS_ZOD)];
%             interf_noise     = 10*log10(sum(10.^(RSRPt/10)) + Noise_mW); 
%             Ue.SINR_geometry = (scenario.Tx_power+pl)-interf_noise;
%             link             = link_list(link_id); 
            [~,~,~,SV(:,n1)] = Ue.get_channel(BSsector_list(link_id),t);
            RSV(n1) = max(SV(:,n1))-min(SV(:,n1));
            Ue.attachBS_pos  = link.BS_pos_wrap;
            BSsector_list(link_id).addUE(Ue);

            waitbar(nn/numel(UE_list),wb,['(',scenario.name,', ',num2str(frequency/1e9),'GHz)  ',num2str(floor(nn/numel(UE_list)*1000)/10),' %'])
        end
        close(wb);
        
%%
        if plot_flag > 1
            for t = 1:numel(BSsector_list)
                BSsector_list(t).plot_links();
            end
        end

        %% result
        clear result;
        result.Pr           = (1:numel(CouplingLoss))/numel(CouplingLoss)*100;
        result.CouplingLoss = sort(CouplingLoss);
        result.SIR          = sort(Geometry);
        result.ZSA_ZSD      = sort(ZSA_ZSD,2);
        result.SV           = sort(sort(SV,1),2);
        result.RSV          = sort(RSV);

        figure(2);
        plot(result.CouplingLoss, result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of CouplingLoss');xlabel('CouplingLoss (dB)');ylabel('CDF'); xlim([-160, -40])

        figure(3);
        plot(result.SIR, result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of geometry SIR');xlabel('geometry SIR (dB)');ylabel('CDF'); xlim([-10, 30])
        
        figure(4);
        plot(result.ZSA_ZSD(1,:), result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of ZSA');xlabel('ZSA (degree)');ylabel('CDF'); xlim([0, 50])

        figure(5);
        plot(result.ZSA_ZSD(2,:), result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of ZSD');xlabel('ZSD (degree)');ylabel('CDF'); xlim([0, 50])

        figure(6);
        plot(result.SV(2,:), result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of 1st angular value');xlabel('1st angular value');ylabel('CDF'); xlim([0, 30])
        
        figure(7);
        plot(result.SV(1,:), result.Pr, 'linewidth',1.5);grid on; hold on;
        title('CDF of 2nd angular value');xlabel('2nd angular value');ylabel('CDF'); xlim([-25, 25])

        figure(8);
        plot(result.RSV, result.Pr, 'linewidth',1.5);grid on; hold on;
        title('CDF of ratio of singular values');xlabel('Ratio of singular values');ylabel('CDF'); xlim([0, 50])
        
        fname = ['.\result\Phase2_36873_result','_',today,'\',scenario.name,'_K=',num2str(BSant_params.panel.Kv),'_M=',num2str(BSant_params.panel.M),'_',num2str(frequency/1e9),'GHz.mat'];
        save(fname,'result');
    end
end

for fig = 2:8
    figure(fig);
    legend('3D-UMa (K=1, M=2)','3D-UMi (K=1, M=2)','3D-UMa (K=M=10)','3D-UMi (K=M=10)','location','southeast');
end