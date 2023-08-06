clear
cd ..
TR = '38901';
date = '2020-11-10';
if ~exist(['.\results\TR',TR,'\InF_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{ 'Geometry_SIR'},{ 'DS'},{ 'ASD'},{ 'ZSD'},{ 'ASA'},{ 'ZSA'}];
sheetName = ["Sub-scenario 1 - 3.5 GHz","Sub-scenario 1 - 28 GHz","Sub-scenario 2 - 3.5 GHz","Sub-scenario 2 - 28 GHz","Sub-scenario 3 - 3.5 GHz","Sub-scenario 3 - 28 GHz","Sub-scenario 4 - 3.5 GHz","Sub-scenario 4 - 28 GHz"];
f = [];
for fign = 1:numel(sheetName)*numel(picname)
    f = [f,figure];
end

%%
freq = [3.5e9, 28e9];
InFcase = {'SL','DL','SH','DH'};
for ii = 1:4
    for jj = 1:2
        frequency = freq(jj);
        fname = ['.\results\TR',TR,'\InF_result_',date,'\InF_',InFcase{ii},'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
        figure(f(((ii-1)*2+jj-1)*8+1));
        plot(result.CouplingLoss, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('CouplingLoss (dB)');ylabel('CDF'); 
        figure(f(((ii-1)*2+jj-1)*8+2));
        plot(result.SINR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
        xlabel('Geometry SINR (dB)');ylabel('CDF'); xlim([floor(min(result.SINR)/10)*10, 30])
        figure(f(((ii-1)*2+jj-1)*8+3));
        plot(result.SIR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF'); xlim([floor(min(result.SIR)/10)*10, 30])
        figure(f(((ii-1)*2+jj-1)*8+4));
        plot(result.DS, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of DS');
        xlabel('DS (ns)');ylabel('CDF');
        figure(f(((ii-1)*2+jj-1)*8+5));
        plot(result.SA4(3,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ASD');
        xlabel('ASD (degree)');ylabel('CDF'); xlim([0, 120])
        figure(f(((ii-1)*2+jj-1)*8+6));
        plot(result.SA4(1,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ZSD');
        xlabel('ZSD (degree)');ylabel('CDF'); xlim([0, 50])
        figure(f(((ii-1)*2+jj-1)*8+7));
        plot(result.SA4(4,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ASA');
        xlabel('ASA (degree)');ylabel('CDF'); xlim([0, 120])
        figure(f(((ii-1)*2+jj-1)*8+8));
        plot(result.SA4(2,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ZSA');
        xlabel('ZSA (degree)');ylabel('CDF'); xlim([0, 50])
    end
end

for k = 1:numel(sheetName)
    data3GPP = data.importfile_InF(".\Docs\R1-1909704_InF_calibration\R1-1909704 Indoor industrial channel model calibration results.xlsx", sheetName(k));
    for kk = 1:numel(picname)
        figure(f((k-1)*numel(picname)+kk));
%         figure(f(kk));
        plot(data3GPP(:,(kk*7)-1),0:100,'r-','linewidth',1.5); hold on;
        plot(data3GPP(:,((kk*7)-5):((kk*7)-2)),0:100,'-','color',[0.7,0.7,0.7],'linewidth',2); hold on;
        plot(data3GPP(:,(kk*7)-1),0:100,'r-','linewidth',1.5); hold on;
    end
end
%%
freq = [3.5e9, 28e9];
InFcase = {'SL','DL','SH','DH'};
for ii = 1:4
    for jj = 1:2
        frequency = freq(jj);
        fname = ['.\results\TR',TR,'\InF_result_',date,'\InF_',InFcase{ii},'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
        figure(f(((ii-1)*2+jj-1)*8+1));
        plot(result.CouplingLoss, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('CouplingLoss (dB)');ylabel('CDF'); 
        figure(f(((ii-1)*2+jj-1)*8+2));
        plot(result.SINR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
        xlabel('Geometry SINR (dB)');ylabel('CDF'); xlim([floor(min(result.SINR)/10)*10, 30])
        figure(f(((ii-1)*2+jj-1)*8+3));
        plot(result.SIR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF'); xlim([floor(min(result.SIR)/10)*10, 30])
        figure(f(((ii-1)*2+jj-1)*8+4));
        plot(result.DS, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of DS');
        xlabel('DS (ns)');ylabel('CDF');
        figure(f(((ii-1)*2+jj-1)*8+5));
        plot(result.SA4(3,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ASD');
        xlabel('ASD (degree)');ylabel('CDF'); xlim([0, 120])
        figure(f(((ii-1)*2+jj-1)*8+6));
        plot(result.SA4(1,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ZSD');
        xlabel('ZSD (degree)');ylabel('CDF'); xlim([0, 50])
        figure(f(((ii-1)*2+jj-1)*8+7));
        plot(result.SA4(4,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ASA');
        xlabel('ASA (degree)');ylabel('CDF'); xlim([0, 120])
        figure(f(((ii-1)*2+jj-1)*8+8));
        plot(result.SA4(2,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ZSA');
        xlabel('ZSA (degree)');ylabel('CDF'); xlim([0, 50])
    end
end
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)*8
    figure(f(a));
    legend('InF','InF, mean of R1-1909704','InF, R1-1909704','location','southeast');
    saveas(f(a),['.\results\TR',TR,'\InF_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_InF_',num2str(a),'.png']);
end