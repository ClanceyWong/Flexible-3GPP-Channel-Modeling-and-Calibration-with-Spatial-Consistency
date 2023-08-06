plot_flag   = 1;
set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type
version = 38901;

N0_dBm      = -174;
N0_linear   = (10^(N0_dBm/10));  % mW/Hz 

sec_per_BS  = 1;              % the number of sector per BS.
ue_per_sec  = 90;             % the number of ue per sector.
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
UEant_params.panel.P            = 1;
UEant_params.panel.X_pol        = [0]; 
UEant_params.panel.ele_downtilt = 90;
UEant_params.panel.ele_panning  = 0;
beta                            = 0;


today = char(datetime('today'));
% today = '2020-08-06';
if ~exist(['.\results\TR38901\InF_result','_',today],'file') 
    mkdir(['.\results\TR38901\InF_result','_',today]);
end 

bandwidth   = 1e8;           % 100 MHz
InFcase     = {'SL','DL','SH','DH'};

for frequency = [3.5e9, 28e9]
    if frequency == 3.5e9
        bandwidth = 1e8;
    else
        bandwidth   = 1e8;           % 100 MHz
    end
    for case_num = 1:4
        BSant_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
        BSant_params.array.Ng           = 1;
        BSant_params.array.dg_H         = 2.5;
        BSant_params.array.dg_V         = 2.5;
        BSant_params.panel.M            = 1;                 %  the number of antenna elements with the same polarization in each column.
        BSant_params.panel.N            = 1;
        BSant_params.panel.Kv           = 1;                  %  K = M or 1
        BSant_params.panel.Kh           = 1;
        BSant_params.panel.d_H          = 0.5; 
        BSant_params.panel.d_V          = 0.5; 
        BSant_params.panel.P            = 1;
        BSant_params.panel.X_pol        = [0]; 
        BSant_params.panel.ele_downtilt = 90;
        BSant_params.panel.ele_panning  = 0;

        scenario = scenarios.InF(InFcase{case_num});
        scenario.BW = bandwidth;

        Noise_mW    = 10^((10*log10(N0_linear * bandwidth) + 9)/10);

        % drop BS
        [BS, BSsector_list]    = network_layout.Drop_BaseStation(scenario, BSant_params, sec_per_BS, frequency, version, plot_flag); % , 'SpatialConsistency'
        % drop UE
        if strcmp(scenario.name,'3D-RMa') ||strcmp(scenario.name,'3D-UMa')||strcmp(scenario.name,'3D-UMi')
            [UE_list, ue_pos_list] = network_layout.Drop_UE(BSsector_list, UEant_params, ue_per_sec, scenario.R, scenario.min_d, beta, plot_flag);
        elseif strcmp(scenario.name,'3D-InH') 
            [UE_list, ue_pos_list] = network_layout.Drop_UE_InH(BSsector_list, UEant_params, ue_per_sec, scenario.min_d, [scenario.x_range;scenario.y_range], beta, plot_flag);
        elseif strcmp(scenario.name,'InF')
            [UE_list, ue_pos_list] = network_layout.Drop_UE_InF(BSsector_list, UEant_params, ue_per_sec, scenario.min_d, [scenario.x_range;scenario.y_range], beta, plot_flag);
        end
        [~, indx]              = sort(ue_pos_list(:,3));

        %% create all links
        port0 = 0;
        CouplingLoss = zeros(1,numel(UE_list)); 
        RSRP         = zeros(1,numel(UE_list)); 
        Geometry     = zeros(1,numel(UE_list));
        GeometryN    = zeros(1,numel(UE_list));
        DS           = zeros(1,numel(UE_list));
        ds           = zeros(1,numel(UE_list));
        SA4          = zeros(4,numel(UE_list));
        as           = zeros(4,numel(UE_list));
        wb           = waitbar(0,'creating all links...');
        for nn = 1:numel(UE_list)
            n1        = indx(nn);
%             if UE_list(n1).h_UT ~= UE_list(indx(max(nn-1,1))).h_UT
%                 for fl = 1:numel(BS)
%                     BS(fl).getRandnGrid(version);
%                 end
%             end
            link_list = [];
            Ue           = UE_list(n1);
            pathloss     = zeros(size(BSsector_list)); 
            couplingloss = zeros(size(BSsector_list));
            RSRPt        = zeros(size(BSsector_list));
            for n2 = 1:numel(BS)
                site      = BS(n2);
                link      = links.Link38901(site, Ue, scenario, true, tab);
                for n3 = 1:sec_per_BS
                    [RSRPt(n2), couplingloss(n2),pathloss(n2)] = Ue.RSRP_calc(BSsector_list(site.sector(n3)),link,port0);
                end
                link_list = [link_list;link]; %#ok<AGROW>
%                 RSRP      = link.RSRP_cal();
%                 PL(n2)    = -link.PL + link.SF + RSRP{1}(1);
%                 link.CouplingLoss = PL(n2);
            end
            [pl, link_id]    = max(RSRPt);
            link             = link_list(ceil(link_id/3));
            BSsector_list(link_id).init_link_params(link);
            RSRP(n1)         = RSRPt(link_id); 
            CouplingLoss(n1) = couplingloss(link_id); 
            RSRPt(link_id)   = -inf;
            interf           = 10*log10(sum(10.^(RSRPt/10)));
            Geometry(n1)     = RSRP(n1)-interf;
            Ue.couplingloss  = CouplingLoss(n1);
            Ue.SIR_geometry  = Geometry(n1);
            interf_noise     = 10*log10(sum(10.^(RSRPt/10)) + Noise_mW); 
            Ue.SINR_geometry = RSRP(n1)-interf_noise;
            GeometryN(n1)    = Ue.SINR_geometry;
            link             = link_list(link_id); 
%             sf = [sf, link.SF];
            Ue.attachBS_pos  = link.BS_pos_wrap;
%             Ue.attachedSector= link.sector;
%             CouplingLoss(n1) = link.CouplingLoss; 
            DS(n1)       = link.calc_delay_spread(link.tau_n_LOS, link.Pn );
            ds(n1)       = link.DS;
            SA4(:,n1)    = [link.calc_angular_spreads_(link.theta_n_m_ZOD,link.Pn,link.theta_LOS_ZOD); ...
                            link.calc_angular_spreads_(link.theta_n_m_ZOA,link.Pn,link.theta_LOS_ZOA); ...
                            link.calc_angular_spreads_(link.phi_n_m_AOD,link.Pn,link.phi_LOS_AOD); ...
                            link.calc_angular_spreads_(link.phi_n_m_AOA,link.Pn,link.phi_LOS_AOA)];
            as(:,n1)     = [link.calc_angular_spreads(link.theta_n_m_ZOD,link.Pn,link.theta_LOS_ZOD); ...
                            link.calc_angular_spreads(link.theta_n_m_ZOA,link.Pn,link.theta_LOS_ZOA); ...
                            link.calc_angular_spreads(link.phi_n_m_AOD,link.Pn,link.phi_LOS_AOD); ...
                            link.calc_angular_spreads(link.phi_n_m_AOA,link.Pn,link.phi_LOS_AOA)];
%             link.sector.addUE(Ue);

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
        result.SINR         = sort(GeometryN);
        result.DS           = sort(DS)*1e9;
        result.ds           = sort(ds)*1e9;
        result.SA4          = sort(SA4,2);
        result.as          = sort(as,2);
%         result.SINR         = sort([UE_list(:).SINR_geometry]);
        
%         figure(100);
%         plot(sort(sf), result.Pr, 'linewidth',1.5);grid on;hold on;

        figure(2);
        plot(result.CouplingLoss, result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of CouplingLoss');xlabel('CouplingLoss (dB)');ylabel('CDF');% xlim([-180, -40])

        figure(3);
        plot(result.SIR, result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of geometry SIR');xlabel('geometry SIR (dB)');ylabel('CDF');% xlim([-30, 30])

        figure(4);
        plot(result.SINR, result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of geometry SINR');xlabel('geometry SINR (dB)');ylabel('CDF');% xlim([-30, 30])
        
        figure(5);
        plot(result.DS, result.Pr, 'linewidth',1.5);grid on;hold on;
        plot(result.ds, result.Pr, '--', 'linewidth',1.5);grid on;hold on;
        title('CDF of Delay spread');xlabel('DS (ns)');ylabel('CDF'); xlim([0, 500])
        
        figure(6);
        plot(result.SA4(1,:), result.Pr, 'linewidth',1.5);grid on; hold on;
        plot(result.as(1,:), result.Pr,'--', 'linewidth',1.5);grid on; hold on;
        title('CDF of ZSD');xlabel('ZSD (degree)');ylabel('CDF');xlim([0, 70])

        figure(7);
        plot(result.SA4(2,:), result.Pr, 'linewidth',1.5);grid on; hold on;
        plot(result.as(2,:), result.Pr,'--', 'linewidth',1.5);grid on; hold on;
        title('CDF of ZSA');xlabel('ZSA (degree)');ylabel('CDF');xlim([0, 70])
        
        figure(8);
        plot(result.SA4(3,:), result.Pr, 'linewidth',1.5);grid on; hold on;
        plot(result.as(3,:), result.Pr,'--', 'linewidth',1.5);grid on; hold on;
        title('CDF of ASD');xlabel('ASD (degree)');ylabel('CDF');xlim([0, 120])

        figure(9);
        plot(result.SA4(4,:), result.Pr, 'linewidth',1.5);grid on; hold on;
        plot(result.as(4,:), result.Pr,'--', 'linewidth',1.5);grid on; hold on;
        title('CDF of ASA');xlabel('ASA (degree)');ylabel('CDF');xlim([0, 120])

        fname = ['.\results\TR38901\InF_result','_',today,'\',scenario.name,'_',InFcase{case_num},'_',num2str(frequency/1e9),'GHz.mat'];
        save(fname,'result');
    end
end

% for fig = 2:8
%     figure(fig);
%     legend('3D-UMa','3D-UMi','3D-InH','location','southeast');
% end