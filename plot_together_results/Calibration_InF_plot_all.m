clear
cd ..
TR = '38901';
date = '2020-11-14';
if ~exist(['.\results\TR',TR,'\InF_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{ 'Geometry_SIR'},{ 'DS'},{ 'ASD'},{ 'ZSD'},{ 'ASA'},{ 'ZSA'}];
sheetName = ["Sub-scenario 1 - 3.5 GHz","Sub-scenario 1 - 28 GHz","Sub-scenario 2 - 3.5 GHz","Sub-scenario 2 - 28 GHz","Sub-scenario 3 - 3.5 GHz","Sub-scenario 3 - 28 GHz","Sub-scenario 4 - 3.5 GHz","Sub-scenario 4 - 28 GHz"];
f = [];
for fign = 1:numel(picname)
    f = [f,figure];
end

%%
hr = [];
for k = 1:numel(sheetName)
    data3GPP = data.importfile_InF(".\Docs\R1-1909704_InF_calibration\R1-1909704 Indoor industrial channel model calibration results.xlsx", sheetName(k));
    for kk = 1:numel(picname)
        figure(f(kk));
%         figure(f(kk));
%         plot(data3GPP(:,(kk*7)-1),0:100,'r-','linewidth',1.5); hold on;
        plot(data3GPP(:,((kk*7)-5):((kk*7)-2)),0:100,'-','color',[0.85,0.85,0.85],'linewidth',3); hold on;
    end
end
 marker = 'sxo+';
 msize = [11,12,12,12];
 color = [0 0.447058823529412 0.741176470588235;...
     0.850980392156863 0.325490196078431 0.098039215686274;...
     0.4 0.6 0.12;...
     0.494117647058824 0.184313725490196 0.556862745098039];
for k = 1:numel(sheetName)
    data3GPP = data.importfile_InF(".\Docs\R1-1909704_InF_calibration\R1-1909704 Indoor industrial channel model calibration results.xlsx", sheetName(k));
    ind = floor((k-1)/2)+1;
    inc = mod(k-1,2)+1;
    for kk = 1:numel(picname)
        figure(f(kk));
%         figure(f(kk));
%         plot(data3GPP(:,(kk*7)-1),0:100,'r-','linewidth',1.5); hold on;
        hrr = plot(data3GPP(:,(kk*7)-1),0:100,'MarkerIndices',[5 15 25 35 45 55 65 75 85 95]+1,...% 'MarkerSize',msize(inc), 'Marker',marker(inc),...
            'LineWidth',1.5,...
            'LineStyle','-','Color',color(ind,:)); hold on;
        hr = [hr,hrr];
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
        [~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
        
        figure(f(1));
        plot(result.CouplingLoss, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 
        figure(f(2));
        plot(result.SINR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
        xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); xlim([floor(min(result.SINR)/10)*10, 30])
        figure(f(3));
        plot(result.SIR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF (%)'); xlim([floor(min(result.SIR)/10)*10, 30])
        figure(f(4));
        plot(result.DS, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of DS');
        xlabel('DS (ns)');ylabel('CDF (%)');
        figure(f(5));
        plot(result.SA4(3,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',2,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of ASD');
        xlabel('ASD (degree)');ylabel('CDF (%)'); xlim([0, 120])
        figure(f(6));
        plot(result.SA4(1,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',2,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of ZSD');
        xlabel('ZSD (degree)');ylabel('CDF (%)'); xlim([0, 50])
        figure(f(7));
        plot(result.SA4(4,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',2,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of ASA');
        xlabel('ASA (degree)');ylabel('CDF (%)'); xlim([0, 120])
        figure(f(8));
        plot(result.SA4(2,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(jj),...
            'Marker',marker(jj),...
            'LineWidth',2,...
            'LineStyle','--','Color',color(ii,:));grid on;hold on;
%         title('CDF of ZSA');
        xlabel('ZSA (degree)');ylabel('CDF (%)'); xlim([0, 50])
    end
end
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)% *8
    figure(f(a));
    m7 = plot(0,-10,'k-',...
    'LineWidth',2); hold on;
m9 = plot(0,-10,'k--',...
    'LineWidth',2);

m1 = plot(0,-10,'-',...
    'Color',color(1,:),...
    'LineWidth',2); hold on;
m2 = plot(0,-10,'-',...
    'Color',color(2,:),...
    'LineWidth',2); hold on;
m3 = plot(0,-10,'-',...
    'Color',color(3,:),...
    'LineWidth',2);
m6 = plot(0,-10,'-',...
    'Color',color(4,:),...
    'LineWidth',2);

m4 = plot([0,0],[-10,-20],'k',...
    'LineWidth',1.5, 'Marker',marker(1)); hold on;
m5 = plot([0,0],[-10,-20],'k',...
    'LineWidth',1.5, 'Marker',marker(2)); hold on;

legend([m7,m9,m2,m3,m6,m1,m4,m5], 'Ref [50]','TR 38.901','DL','SH','DH','SL','3.5 GHz','28 GHz','fontsize',14)
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
    ylim([0 100]);
end