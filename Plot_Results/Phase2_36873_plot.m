clear
cd ..
date = '2020-07-23';
if ~exist(['.\result\Phase2_36873_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

% index: 1:Geometry SINR 
%        2:ZSD
%        3:ZSA
%        4:Ratio singular value
%        5:Largest singular value
%        6:Smallest singular value
%        7:Coupling loss
index = [7,1,2,3,4,5,6];
f1 = figure; f2 = figure; f3 = figure; f4 = figure; f5 = figure; f6 = figure; f7 = figure; f = [f1,f2,f3,f4,f5,f6,f7];
picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{'ZSD'},{'ZSA'},{'Ratio_singular_value'},{'Largest_singular_value'},{'Smallest_singular_value'}];

sheetName = ["3D-UMa (K=1,M=2)","3D-UMa (K=M=10)","3D-UMi (K=1,M=2)","3D-UMi (K=M=10)"];
for k = 1:4
    data3GPP = data.importfile_p2(".\data_3GPP\NOKNET_3DPhase2CalibrationGeogrDistance_v59.xls", sheetName(k), [29, 129]);
    for kk = 1:numel(picname)
        figure(f(kk));
        plot(data3GPP(:,index(kk)),0:100,'--','linewidth',1.5); hold on;
    end
end
freq = 2;
for ii = 1:4
    if ii == 1
        name = '3D-UMa';    
        M = 2;  K = 1;   
    elseif ii == 2
        name = '3D-UMa';    
        M = 10; K = M;     
    elseif ii == 3
        name = '3D-UMi';  
        M = 2;  K = 1;    
    elseif ii == 4
        name = '3D-UMi';      
        M = 10; K = M;    
    end
    fname = ['.\result\Phase2_36873_result_',date,'\',name,'_K=',num2str(K),'_M=',num2str(M),'_',num2str(freq),'GHz.mat'];
    load(fname,'result');
    figure(f1);
    plot(result.CouplingLoss, result.Pr, 'linewidth',2);grid on;hold on;
%     title('CDF of CouplingLoss');
    xlabel('CouplingLoss (dB)');ylabel('CDF');xlim([-160, -40])
    figure(f2);
    plot(result.SIR, result.Pr, 'linewidth',2);grid on;hold on;
%     title('CDF of Geometry');
    xlabel('Geometry (dB)');ylabel('CDF');xlim([-10, 30])
    figure(f3);
    plot(result.ZSA_ZSD(2,:), result.Pr, 'linewidth',1.5);grid on; hold on;
%     title('CDF of ZSD');
    xlabel('ZSD (degree)');
    ylabel('CDF');xlim([0, 50])
    figure(f4);
    plot(result.ZSA_ZSD(1,:), result.Pr, 'linewidth',1.5);grid on; hold on;
%     title('CDF of ZSA');
    xlabel('ZSA (degree)');ylabel('CDF');xlim([0, 50])
    figure(f5);
    plot(result.RSV, result.Pr, 'linewidth',1.5);grid on; hold on;
%     title('CDF of Ratio of singular values');
    xlabel('Ratio of singular values (dB)');ylabel('CDF');xlim([0, 50]);
    figure(f6);
    plot(result.SV(2,:), result.Pr, 'linewidth',1.5);grid on; hold on;
%     title('CDF of Largest singular value');
    xlabel('Largest singular value (dB)');ylabel('CDF');xlim([-5,30]);
    figure(f7);
    plot(result.SV(1,:), result.Pr, 'linewidth',1.5);grid on; hold on;
%     title('CDF of Smallest singular value');
    xlabel('Smallest singular value (dB)');ylabel('CDF');xlim([-25,25]);
end

for a = 1:numel(picname)
    figure(f(a));
    legend('3D-UMa (K=1,M=2), R1-143469','3D-UMa (K=M=10), R1-143469','3D-UMi (K=1,M=2), R1-143469','3D-UMi (K=M=10), R1-143469','3D-UMa (K=1,M=2)','3D-UMa (K=M=10)','3D-UMi (K=1,M=2)','3D-UMi (K=M=10)','location','southeast');
    saveas(f(a),['.\result\Phase2_36873_result_',date,'\',cell2mat(picname(a)),'.png']);
end