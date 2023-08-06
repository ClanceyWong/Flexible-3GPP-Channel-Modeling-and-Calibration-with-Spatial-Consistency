clear
cd ..
TR = '38900';
date = '2020-10-19';
config = '1';
if ~exist(['.\results\TR',TR,'\spatial_consistency_result_config',config,'_3_6_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

if strcmp(config,'1')
    picname = [{'delay'},{ 'AOA'},{'LOS state'},{'CIR'}];
    sheetName = ["Config1 (metric 3-6)"];
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
fname = ['.\results\TR',TR,'\spatial_consistency_result_config',config,'_3_6_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
load(fname,'result');
figure(f(1));
plot(result.distance, result.xcorr(1,:)*100,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('distane');ylabel('xcorr-delay'); 

figure(f(2));
%         title('CDF of Geometry SIR');
plot(result.distance, result.xcorr(2,:)*100,'b-', 'linewidth',2);grid on;hold on;
xlabel('distane');ylabel('xcorr-AOA'); 

figure(f(3));
%         title('CDF of Geometry SIR');
plot(result.distance, result.xcorr(3,:)*100,'b-', 'linewidth',2);grid on;hold on;
xlabel('distane');ylabel('xcorr-LOS state'); 
% 
figure(f(4));
% %         title('CDF of Geometry CIR');
plot(result.distance, result.xcorr(4,:),'b-', 'linewidth',2);grid on;hold on;
xlabel('distane');ylabel('xcorr-CIR'); 

data3GPP = data.importfile(".\Docs\R1-1609785_addition_calibration\Phase3SpatialConsistency_v9_Xinwei.xlsx", sheetName(1),[29, 159]);
for kk = 1:numel(picname)
    figure(f(kk));
%         figure(f(kk));
    plot(0:130,data3GPP(:,(kk*20)),'r-','linewidth',1.5); hold on;
    plot(0:130,data3GPP(:,((kk-1)*20+1):(kk*20-1)),'-','color',[0.7,0.7,0.7],'linewidth',2); hold on;
    plot(0:130,data3GPP(:,(kk*20)),'r-','linewidth',1.5); hold on;
end
%%
figure(f(1));
plot(result.distance, result.xcorr(1,:)*100,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('distane');ylabel('xcorr-delay'); 

figure(f(2));
%         title('CDF of Geometry SIR');
plot(result.distance, result.xcorr(2,:)*100,'b-', 'linewidth',2);grid on;hold on;
xlabel('distane');ylabel('xcorr-AOA'); 

figure(f(3));
%         title('CDF of Geometry SIR');
plot(result.distance, result.xcorr(3,:)*100,'b-', 'linewidth',2);grid on;hold on;
xlabel('distane');ylabel('xcorr-LOS state'); 

figure(f(4));
% %         title('CDF of Geometry CIR');
plot(result.distance, result.xcorr(4,:),'b-', 'linewidth',2);grid on;hold on;
xlabel('distane');ylabel('xcorr-CIR'); 
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    legend('UMi','UMi, mean of R1-1609785','UMi, R1-1609785','location','northeast');
    saveas(f(a),['.\results\TR',TR,'\spatial_consistency_result_config',config,'_3_6_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
 end