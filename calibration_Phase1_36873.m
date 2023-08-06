plot_flag   = 1;
set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

N0_dBm      = -174;
N0_linear   = (10^(N0_dBm/10));  % mW/Hz 

layer_num   = 3;              % the number of layer of cells.
sec_per_BS  = 3;              % the number of sector per BS.
ue_per_sec  = 50;             % the number of ue per sector.
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
UEant_params.panel.X_pol        = 0; 
UEant_params.panel.ele_downtilt = 90;
UEant_params.panel.ele_panning  = 0;
beta                            = 0;


today = char(datetime('today'));
% today = '2020-07-20';
if ~exist(['.\result\Phase1_36873_result','_',today],'file') 
    mkdir(['.\result\Phase1_36873_result','_',today]);
end 

bandwidth   = 1e7;           % 10 MHz
for frequency = 2e9 %
    for case_num = 1:4
        if case_num == 1 || case_num == 2
            BSant_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
            BSant_params.array.Ng           = 1;
            BSant_params.array.dg_H         = 2.5;
            BSant_params.array.dg_V         = 2.5;
            BSant_params.panel.M            = 10;                 %  the number of antenna elements with the same polarization in each column.
            BSant_params.panel.Kv           = 10;                  %  K = M or 1
            BSant_params.panel.Kh           = 1;
            BSant_params.panel.N            = 1;
            BSant_params.panel.d_H          = 0.5; 
            BSant_params.panel.d_V          = 0.5; 
            BSant_params.panel.P            = 1;
            BSant_params.panel.X_pol        = 0; 
            BSant_params.panel.ele_downtilt = 102;
            BSant_params.panel.ele_panning  = 0;
%             BSant_params.panel.pol_model    = 'model-1';
        else
            BSant_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
            BSant_params.array.Ng           = 1;
            BSant_params.array.dg_H         = 2.5;
            BSant_params.array.dg_V         = 2.5;
            BSant_params.panel.M            = 1;                 %  the number of antenna elements with the same polarization in each column.
            BSant_params.panel.Kv           = 1;                  %  K = M or 1
            BSant_params.panel.Kh           = 1;
            BSant_params.panel.N            = 1;
            BSant_params.panel.d_H          = 0.5; 
            BSant_params.panel.d_V          = 0.5; 
            BSant_params.panel.P            = 1;
            BSant_params.panel.X_pol        = 0; 
            BSant_params.panel.ele_downtilt = 102;
            BSant_params.panel.ele_panning  = 0;
        end

        if case_num == 1 || case_num == 3
            scenario = scenarios.UMa_3D(layer_num);
            scenario.BW = bandwidth;
        elseif case_num == 2  || case_num == 4
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
        CouplingLoss  = zeros(1,numel(UE_list)); 
        Geometry_SIR  = zeros(1,numel(UE_list));
        Geometry_SINR = zeros(1,numel(UE_list));
        LOS_ZOD       = zeros(1,numel(UE_list));
        wb            = waitbar(0,'creating all links...');
        for nn = 1:numel(UE_list)
            n1        = indx(nn);
%             if UE_list(n1).h_UT ~= UE_list(indx(max(nn-1,1))).h_UT
%                 for fl = 1:numel(BS)
%                     BS(fl).getRandnGrid();
%                 end
%             end
            link_list = [];
            Ue        = UE_list(n1);
            PL        = zeros(size(BSsector_list));
            couplingloss = zeros(size(BSsector_list));
            RSRPt        = zeros(size(BSsector_list));
            for n2 = 1:numel(BS)
                site      = BS(n2);
                link      = links.Link(site, Ue, scenario, true, tab);
                link_list = [link_list;link]; %#ok<AGROW>
                for n3 = 1:sec_per_BS
                    array_gain = BSsector_list(site.sector(n3)).antenna.array_gain(link);
                    couplingloss(3*(n2-1)+n3) = -link.PL + link.SF + array_gain{1}(1);
                    RSRPt(3*(n2-1)+n3)        = -link.PL + link.SF + array_gain{1}(1) + site.Tx_power;
                end
            end
            [pl, link_id]    = max(couplingloss);
            link             = link_list(ceil(link_id/3));
            BSsector_list(link_id).init_link_params(link);
            RSRP             = RSRPt(link_id); 
            RSRPt(link_id)   = -inf;
            interf           = 10*log10(sum(10.^(RSRPt/10)));
            interf_noise     = 10*log10(sum(10.^(RSRPt/10)) + Noise_mW); 
            LOS_ZOD(n1)      = link.theta_LOS_ZOD;
            CouplingLoss(n1) = couplingloss(link_id); 
            Geometry_SIR(n1) = RSRP - interf;
            Geometry_SINR(n1)= RSRP - interf_noise;
            Ue.couplingloss  = CouplingLoss(n1);
            Ue.SIR_geometry  = Geometry_SIR(n1);
            Ue.SINR_geometry = Geometry_SINR(n1);
            
            BSsector_list(link_id).addUE(Ue);
            waitbar(nn/numel(UE_list),wb,['(',scenario.name,', ',num2str(frequency/1e9),'GHz)  ',num2str(floor(nn/numel(UE_list)*100)),' %'])
        end
        close(wb);

        if plot_flag > 1
            for t = 1:numel(BSsector_list)
                BSsector_list(t).plot_links();
            end
        end

        %% result
        clear result;
        result.Pr           = (1:numel(CouplingLoss))/numel(CouplingLoss)*100;
        result.CouplingLoss = sort(CouplingLoss);
        result.SIR          = sort(Geometry_SIR);
        result.SINR         = sort(Geometry_SINR);
        result.LOS_ZOD      = sort(LOS_ZOD);
        
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
        plot(result.LOS_ZOD, result.Pr, 'linewidth',1.5);grid on;hold on;
        title('CDF of LOS ZOD');xlabel('LOS ZOD (degree)');ylabel('CDF');% xlim([-30, 30])
        

        fname = ['.\result\Phase1_36873_result','_',today,'\',scenario.name,'_K=M=',num2str(BSant_params.panel.Kv),'_',num2str(frequency/1e9),'GHz.mat'];
        save(fname,'result');
    end
end

for fig = 2:5
    figure(fig);
    legend('3D-UMa (K=M=10)','3D-UMi (K=M=10)','3D-UMa (K=M=1)','3D-UMi (K=M=1)','location','southeast');
end