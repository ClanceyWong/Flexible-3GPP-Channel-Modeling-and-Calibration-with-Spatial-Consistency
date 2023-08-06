clear
cd ..
TR   = '38900';
date = '2020-08-18';
if ~exist(['.\results\TR',TR,'\large_scale_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type
f1 = figure; f2 = figure; f3 = figure; 
f4 = figure; f5 = figure; f6 = figure;
f7 = figure; f8 = figure; f9 = figure;
f = [f1,f2,f3,f4,f5,f6,f7,f8,f9];
picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{ 'Geometry_SIR'}];

sheetName = ["UMa-6GHz","UMa-30GHz","UMa-70GHz","UMi-6GHz","UMi-30GHz","UMi-70GHz","indoor-6GHz","indoor-30GHz","indoor-70GHz"];
for k = 1:9
    data3GPP = data.importfile_p1(".\Docs\R1-165974_large_scale_calibration\Phase1Calibration_v42_CMCC.xlsx", sheetName(k), [29, 129]);
    for kk = 1:numel(picname)
        figure(f((ceil(k/3)-1)*3+kk));
        plot(data3GPP(:,kk),0:100,'--','linewidth',1.5); hold on;
    end
end
%%
for frequency = [6e9, 3e10,7e10]
    for ii = 1:3
        if ii == 1
            name = '3D-UMa';
        elseif ii == 2
            name = '3D-UMi'; 
        elseif ii == 3
            name = '3D-InH'; 
        end
        fname = ['.\results\TR38901\large_scale_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
        figure(f((ii-1)*3+1));
        plot(result.CouplingLoss, result.Pr, 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('CouplingLoss (dB)');ylabel('CDF'); xlim([floor(min(result.CouplingLoss)/10)*10, -30])
        figure(f((ii-1)*3+2));
        plot(result.SINR, result.Pr, 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SINR');
        xlabel('Geometry SINR (dB)');ylabel('CDF'); xlim([floor(min(result.SINR)/10)*10, 30])
        figure(f((ii-1)*3+3));
        plot(result.SIR, result.Pr,'-', 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF'); xlim([floor(min(result.SIR)/10)*10, 30])
    end
end
%%
loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    legend('UMa-  6 GHz, R1-165974','UMa-30 GHz, R1-165974','UMa-70 GHz, R1-165974','UMa-  6 GHz','UMa-30 GHz','UMa-70 GHz','location',loc{a});
    saveas(f(a),['.\results\TR',TR,'\large_scale_result_',date,'\',cell2mat(picname(a)),'_UMa.png']);
    
    figure(f(a+3));
    legend('UMi-  6 GHz, R1-165974','UMi-30 GHz, R1-165974','UMi-70 GHz, R1-165974','UMi-  6 GHz','UMi-30 GHz','UMi-70 GHz','location',loc{a});
    saveas(f(a+3),['.\results\TR',TR,'\large_scale_result_',date,'\',cell2mat(picname(a)),'_UMi.png']);
    
    figure(f(a+6));
    legend('InH-  6 GHz, R1-165974','InH-30 GHz, R1-165974','InH-70 GHz, R1-165974','InH-  6 GHz','InH-30 GHz','InH-70 GHz','location','southeast');
    saveas(f(a+6),['.\results\TR',TR,'\large_scale_result_',date,'\',cell2mat(picname(a)),'_InH.png']);
end