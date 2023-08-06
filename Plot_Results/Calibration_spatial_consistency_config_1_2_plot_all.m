clear
cd ..
TR = '38900';
date = '2020-10-18';
config = '1';
if ~exist(['.\results\TR',TR,'\spatial_consistency_result_config',config,'_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

if strcmp(config,'1')
    picname = [{'CouplingLoss'},{ 'Geometry_SINR'}];
    sheetName = ["Config1 (metric 1-2)"];
else
    picname = [{'CouplingLoss'},{ 'Geometry_SINR'}];
    sheetName = ["Config2-ProcA"];
end

f = [];
for fign = 1:numel(sheetName)*numel(picname)
    f = [f,figure];
end

%%
frequency = 3e10;
name = '3D-UMi';
fname = ['.\results\TR',TR,'\spatial_consistency_result_config',config,'_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
load(fname,'result');
figure(f(1));
plot(result.CouplingLoss, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('CouplingLoss (dB)');ylabel('CDF'); 

figure(f(2));
%         title('CDF of Geometry SIR');
plot(result.SINR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
xlabel('Geometry SINR (dB)');ylabel('CDF'); 


data3GPP = data.importfile(".\Docs\R1-1609785_addition_calibration\Phase3SpatialConsistency_v9_Xinwei.xlsx", sheetName(1));
for kk = 1:numel(picname)
    figure(f(kk));
%         figure(f(kk));
    plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',1.5); hold on;
    plot(data3GPP(:,((kk-1)*20+1):(kk*20-1)),0:100,'-','color',[0.7,0.7,0.7],'linewidth',2); hold on;
    plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',1.5); hold on;
end
%%
figure(f(1));
plot(result.CouplingLoss, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('CouplingLoss (dB)');ylabel('CDF'); 

figure(f(2));
plot(result.SINR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
xlabel('Geometry SINR (dB)');ylabel('CDF'); 


%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    legend('UMi','UMi, mean of R1-1609785','UMi, R1-1609785','location','southeast');
    saveas(f(a),['.\results\TR',TR,'\spatial_consistency_result_config',config,'_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
 end