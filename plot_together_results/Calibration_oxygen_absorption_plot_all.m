clear
cd ..
TR = '38900';
date = '2020-08-18';
if ~exist(['.\results\TR',TR,'\oxygen_absorption_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{ 'DS'}];
sheetName = ["UMi-60GHz"];
f = [];
for fign = 1:numel(sheetName)*numel(picname)
    f = [f,figure];
end

%%
hr =[];
data3GPP = data.importfile(".\Docs\R1-1609785_addition_calibration\Phase3OxygenAbsorption_v9_Samsung.xlsx", sheetName(1));
for kk = 1:numel(picname)
    figure(f(kk));
%         figure(f(kk));
   % plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',1.5); hold on;
    plot(data3GPP(:,((kk-1)*20+1):(kk*20-1)),0:100,'-','color',[0.85,0.85,0.85],'linewidth',4); hold on;
    hrr = plot(data3GPP(:,(kk*20)),0:100,'r--','linewidth',2); hold on;
    hr = [hr,hrr];
end
%%
frequency = 6e10;
name = '3D-UMi';
fname = ['.\results\TR',TR,'\oxygen_absorption_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
load(fname,'result');
figure(f(1));
h1 = plot(result.CouplingLoss, result.Pr,'b-.', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 

figure(f(2));
h2 = plot(result.SINR, result.Pr,'b-.', 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); 
figure(f(3));
h3 = plot(result.DS, result.Pr,'b-.', 'linewidth',2);grid on;hold on;
%         title('CDF of DS');
xlabel('DS (ns)');ylabel('CDF (%)');
hh0 = [h1,h2,h3];
%%
frequency = 6e10;
name = '3D-UMi';
TR = '38901';
fname = ['.\results\TR',TR,'\oxygen_absorption_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
load(fname,'result');
figure(f(1));
h4 = plot(result.CouplingLoss, result.Pr,'g-', 'linewidth',2); hold on;
%         title('CDF of CouplingLoss');
xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 

figure(f(2));
h5 = plot(result.SINR, result.Pr,'g-', 'linewidth',2); hold on;
%         title('CDF of Geometry SIR');
xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); 
figure(f(3));
h6 = plot(result.DS, result.Pr,'g-', 'linewidth',2); hold on;
%         title('CDF of DS');
xlabel('DS (ns)');ylabel('CDF (%)');
hh1 = [h4,h5,h6];
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    legend([hr(a),hh0(a),hh1(a)],'Ref [49]','TR 38.900','TR 38.901','location','southeast','fontsize',14);
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
%     saveas(f(a),['.\results\TR',TR,'\oxygen_absorption_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
 end