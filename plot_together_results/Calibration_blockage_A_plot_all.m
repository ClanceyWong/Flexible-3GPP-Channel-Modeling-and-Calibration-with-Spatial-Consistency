clear
cd ..
TR = '38900';
date = '2020-08-29';
if ~exist(['.\results\TR',TR,'\blockage_result_A_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{'ASA'}];
sheetName = ["Model A UMi-30GHz"];

f = [];
for fign = 1:numel(sheetName)*numel(picname)
    f = [f,figure];
end

%%
 color = [0 0.447058823529412 0.741176470588235;...
     0.850980392156863 0.325490196078431 0.098039215686274;...
     0.4 0.6 0.12;...
     0.494117647058824 0.184313725490196 0.556862745098039];
hr = [];
data3GPP = data.importfile(".\Docs\R1-1609785_addition_calibration\Phase3Blockage_v11_Samsung.xlsx", sheetName(1));
for kk = 1:numel(picname)
    figure(f(kk));
%         figure(f(kk));
%     plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',1.5); hold on;
    plot(data3GPP(:,((kk-1)*20+1):(kk*20-1)),0:100,'-','color',[0.85,0.85,0.85],'linewidth',3); hold on;
    hrr = plot(data3GPP(:,(kk*20)),0:100,'-','linewidth',2,'color',color(1,:)); hold on;
    hr = [hr,hrr];
end
%%
frequency = 3e10;
name = '3D-UMi';
fname = ['.\results\TR',TR,'\blockage_result_A_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
load(fname,'result');
[~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
figure(f(1));
h1 = plot(result.CouplingLoss, result.Pr,'-.','Marker','s','MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',12, 'linewidth',2,'color',color(2,:));grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 

figure(f(2));
%         title('CDF of Geometry SIR');
h2 = plot(result.SINR, result.Pr,'-.','Marker','s','MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',12, 'linewidth',2,'color',color(2,:));grid on;hold on;
xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); 

figure(f(3));
h3 = plot(result.ASA, result.Pr,'-.','Marker','s','MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',12, 'linewidth',2,'color',color(2,:));grid on;hold on;
xlabel('ASA (degree)');ylabel('CDF (%)'); 
hh0 = [h1,h2,h3];
%%
frequency = 3e10;
name = '3D-UMi';
TR = '38901';
fname = ['.\results\TR',TR,'\blockage_result_A_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
load(fname,'result');
[~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
figure(f(1));
h1 = plot(result.CouplingLoss, result.Pr,'--','Marker','x','MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',12, 'linewidth',2,'color',color(3,:));grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 

figure(f(2));
%         title('CDF of Geometry SIR');
h2 = plot(result.SINR, result.Pr,'--','Marker','x','MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',12, 'linewidth',2,'color',color(3,:));grid on;hold on;
xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); 

figure(f(3));
h3 = plot(result.ASA, result.Pr,'--','Marker','x','MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',12, 'linewidth',2,'color',color(3,:));grid on;hold on;
xlabel('ASA (degree)');ylabel('CDF (%)'); 
hh1 = [h1,h2,h3];
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    legend([hr(a),hh0(a),hh1(a)],'Ref [49]','TR 38.900','TR 38.901','location','southeast','fontsize',14);
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
%     saveas(f(a),['.\results\TR',TR,'\blockage_result_A_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
 end