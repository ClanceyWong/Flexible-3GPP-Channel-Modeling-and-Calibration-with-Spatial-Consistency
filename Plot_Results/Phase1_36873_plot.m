clear
cd ..
date = '2020-07-22';
if ~exist(['.\result\Phase1_36873_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type
f1 = figure; f2 = figure; f3 = figure; f = [f1,f2,f3];
picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{ 'LOS_ZOD'}];

sheetName = ["3D-UMa (K=M=1)","3D-UMa (K=M=10)","3D-UMi (K=M=1)","3D-UMi (K=M=10)"];
for k = 1:4
    data3GPP = data.importfile_p1(".\data_3GPP\3DPhase1CalibrationGeogrDistance_v20.xlsx", sheetName(k), [29, 129]);
    for kk = 1:numel(picname)
        figure(f(kk));
        plot(data3GPP(:,kk),0:100,'--','linewidth',1.5); hold on;
    end
end
freq = 2;
for ii = 1:4
    if ii == 1
        name = '3D-UMa';
        M = 1; 
    elseif ii == 2
        name = '3D-UMa'; 
        M = 10; 
    elseif ii == 3
        name = '3D-UMi'; 
        M = 1; 
    elseif ii == 4
        name = '3D-UMi'; 
        M = 10; 
    end
    fname = ['.\result\Phase1_36873_result_',date,'\',name,'_K=M=',num2str(M),'_',num2str(freq),'GHz.mat'];
    load(fname,'result');
    figure(f1);
    plot(result.CouplingLoss, result.Pr, 'linewidth',2);grid on;hold on;
%     title('CDF of CouplingLoss');
    xlabel('CouplingLoss (dB)');ylabel('CDF');xlim([-160, -40])
    figure(f2);
    plot(result.SINR, result.Pr, 'linewidth',2);grid on;hold on;
%     title('CDF of Geometry SINR');
    xlabel('Geometry SINR (dB)');ylabel('CDF');xlim([-10, 30])
    figure(f3);
    plot(result.LOS_ZOD, result.Pr,'-', 'linewidth',2);grid on;hold on;
%     title('CDF of LOS ZOD');
    xlabel('LOS ZOD (degrees)');ylabel('CDF');xlim([70, 120])
end

for a = 1:numel(picname)
    figure(f(a));
    legend('3D-UMa (K=M=1), R1-143469','3D-UMa (K=M=10), R1-143469','3D-UMi (K=M=1), R1-143469','3D-UMi (K=M=10), R1-143469','3D-UMa (K=M=1)','3D-UMa (K=M=10)','3D-UMi (K=M=1)','3D-UMi (K=M=10)','location','southeast');
    saveas(f(a),['.\result\Phase1_36873_result_',date,'\',cell2mat(picname(a)),'.png']);
end