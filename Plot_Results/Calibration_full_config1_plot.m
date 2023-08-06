clear
cd ..
TR = '38900';
date = '2020-08-03';
if ~exist(['.\results\TR',TR,'\full_config1_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type
f1 = figure; f2 = figure; f3 = figure; f4 = figure; f5 = figure; f6 = figure; f7 = figure;
f8 = figure; f9 = figure; f10 = figure; f11 = figure; f12 = figure; f13 = figure; f14 = figure;
f15 = figure; f16 = figure; f17 = figure; f18 = figure; f19 = figure; f20 = figure; f21 = figure;
f = [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21];
picname = [{'CouplingLoss'},{ 'Geometry_SIR'},{ 'DS'},{ 'ASD'},{ 'ZSD'},{ 'ASA'},{ 'ZSA'}];

sheetName = ["UMa-6GHz","UMa-30GHz","UMa-60GHz","UMa-70GHz","UMi-6GHz","UMi-30GHz","UMi-60GHz","UMi-70GHz","InH-6GHz","InH-30GHz","InH-60GHz","InH-70GHz"];
for k = 1:12
    data3GPP = data.importfile_full(".\Docs\R1-165975_full_calibraton\Phase2Config1Calibration_v35_CMCC.xlsx", sheetName(k), [29, 129]);
    for kk = 1:numel(picname)
        figure(f((ceil(k/4)-1)*7+kk));
%         figure(f(kk));
        plot(data3GPP(:,kk),0:100,'--','linewidth',1.5); hold on;
    end
end
%%
for frequency = [6e9, 3e10, 6e10, 7e10]
    for ii = 1:3
        if ii == 1
            name = '3D-UMa'; xmax = 2000; xm = -40;
        elseif ii == 2
            name = '3D-UMi'; xmax = 2000; xm = -40;
        elseif ii == 3
            name = '3D-InH'; xmax = 200;  xm = -30;
        end
        fname = ['.\results\TR',TR,'\full_config1_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
        figure(f((ii-1)*7+1));
        plot(result.CouplingLoss, result.Pr, 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('CouplingLoss (dB)');ylabel('CDF'); 
        if xm == -30
            xlim([-100, xm])
        else
            xlim([floor(min(result.CouplingLoss)/10)*10, xm])
        end
        figure(f((ii-1)*7+2));
        plot(result.SIR, result.Pr, 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF'); xlim([floor(min(result.SIR)/10)*10, 30])
        figure(f((ii-1)*7+3));
        plot(result.DS, result.Pr,'-', 'linewidth',2);grid on;hold on;
%         title('CDF of DS');
        xlabel('DS (ns)');ylabel('CDF'); xlim([0, xmax])
        figure(f((ii-1)*7+4));
        plot(result.SA4(3,:), result.Pr,'-', 'linewidth',2);grid on;hold on;
%         title('CDF of ASD');
        xlabel('ASD (degree)');ylabel('CDF'); xlim([0, 120])
        figure(f((ii-1)*7+5));
        plot(result.SA4(1,:), result.Pr,'-', 'linewidth',2);grid on;hold on;
%         title('CDF of ZSD');
        xlabel('ZSD (degree)');ylabel('CDF'); xlim([0, 50])
        figure(f((ii-1)*7+6));
        plot(result.SA4(4,:), result.Pr,'-', 'linewidth',2);grid on;hold on;
%         title('CDF of ASA');
        xlabel('ASA (degree)');ylabel('CDF'); xlim([0, 120])
        figure(f((ii-1)*7+7));
        plot(result.SA4(2,:), result.Pr,'-', 'linewidth',2);grid on;hold on;
%         title('CDF of ZSA');
        xlabel('ZSA (degree)');ylabel('CDF'); xlim([0, 50])
    end
end
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    legend('UMa-  6 GHz, R1-165975','UMa-30 GHz, R1-165975','UMa-60 GHz, R1-165975','UMa-70 GHz, R1-165975','UMa-  6 GHz','UMa-30 GHz','UMa-60 GHz','UMa-70 GHz','location','southeast');
    saveas(f(a),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(a)),'_UMa.png']);
    
    figure(f(a+7));
    legend('UMi-  6 GHz, R1-165975','UMi-30 GHz, R1-165975','UMi-60 GHz, R1-165975','UMi-70 GHz, R1-165975','UMi-  6 GHz','UMi-30 GHz','UMi-60 GHz','UMi-70 GHz','location','southeast');
    saveas(f(a+7),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(a)),'_UMi.png']);
    
    figure(f(a+14));
    legend('InH-  6 GHz, R1-165975','InH-30 GHz, R1-165975','InH-60 GHz, R1-165975','InH-70 GHz, R1-165975','InH-  6 GHz','InH-30 GHz','InH-60 GHz','InH-70 GHz','location','southeast');
    saveas(f(a+14),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(a)),'_InH.png']);
end