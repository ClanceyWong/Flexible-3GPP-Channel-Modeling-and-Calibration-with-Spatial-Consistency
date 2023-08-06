clear
cd ..
TR = '38900';
date = '2020-11-10';% 08
if ~exist(['.\results\TR',TR,'\full_config1_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type
set(0,'defaultAxesFontSize',13)                                    % Default Font Type
set(0,'defaultTextFontSize',13)                                    % Default Font Type
 hr = [];
 marker = 'osx+^d';
 linesty = '- --: -.- --: -.- --: -.';
 size = [11,11,11,11];
 color = [0 0.447058823529412 0.741176470588235;...
     0.850980392156863 0.325490196078431 0.098039215686274;...
     0.4 0.6 0.12;...
     0.494117647058824 0.184313725490196 0.556862745098039];

picname = [{'CouplingLoss'},{ 'Geometry_SIR'},{ 'DS'},{ 'ASD'},{ 'ZSD'},{ 'ASA'},{ 'ZSA'}];
sheetName = ["UMa-6GHz","UMa-30GHz","UMa-60GHz","UMa-70GHz","UMi-6GHz","UMi-30GHz","UMi-60GHz","UMi-70GHz","InH-6GHz","InH-30GHz","InH-60GHz","InH-70GHz"];
f = [];
for fign = 1:numel(picname) % numel(sheetName)*
    f = [f,figure];
end
hr = [];
for k = 9:12
    data3GPP = data.importfile(".\Docs\R1-165975_full_calibraton\Phase2Config1Calibration_v35_CMCC.xlsx", sheetName(k));
    for kk = 1:numel(picname)
%         figure(f((k-1)*numel(picname)+kk));
        figure(f(kk));
%         plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',1.5); hold on;
        plot(data3GPP(:,((kk-1)*20+1):(kk*20-1)),0:100,'-','color',[0.85,0.85,0.85],'linewidth',3); hold on;
%         plot(data3GPP(:,(kk*20)),0:100,'-','linewidth',2); hold on;
    end
end
for k = 9:12
    data3GPP = data.importfile(".\Docs\R1-165975_full_calibraton\Phase2Config1Calibration_v35_CMCC.xlsx", sheetName(k));
    ind = floor((k-1)/4)+1;
    inc = mod(k-1,4)+1;
    for kk = 1:numel(picname)
%         figure(f((k-1)*numel(picname)+kk));
        figure(f(kk));
%         plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',1.5); hold on;
%         plot(data3GPP(:,((kk-1)*20+1):(kk*20-1)),0:100,'-','color',[0.85,0.85,0.85],'linewidth',3); hold on;
        hrr = plot(data3GPP(:,(kk*20)),0:100,'MarkerIndices',[5 15 25 35 45 55 65 75 85 95]+1,...%'MarkerSize',size(inc),  'Marker',marker(inc),...
            'LineWidth',2,...
            'LineStyle',linesty((k-1)*2+1:2*k),'Color',color(inc,:)); hold on;
        hr = [hr, hrr];
    end
end
%%
hs = [];
freq = [6e9, 3e10, 6e10, 7e10];
for ii = 3
    for jj = 1:4
        frequency = freq(jj);
        if ii == 1
            name = '3D-UMa'; xmax = 2000; xm = -40;
        elseif ii == 2
            name = '3D-UMi'; xmax = 2000; xm = -40;
        elseif ii == 3
            name = '3D-InH'; xmax = 200;  xm = -30;
        end
        fname = ['.\results\TR',TR,'\full_config1_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
%         figure(f(((ii-1)*4+jj-1)*7+1));
        [~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
        figure(f(1));
        hs1 = plot(result.CouplingLoss, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 
%         if xm == -30
%             xlim([-100, xm])
%         else
%             xlim([floor(min(result.CouplingLoss)/10)*10, xm])
%         end
%         figure(f(((ii-1)*4+jj-1)*7+2));
        figure(f(2));
        hs2 = plot(result.SIR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF (%)'); % xlim([floor(min(result.SIR)/10)*10, 30])
%         figure(f(((ii-1)*4+jj-1)*7+3));
        figure(f((3)));
        hs3 = plot(result.DS, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of DS');
        xlabel('DS (ns)');ylabel('CDF (%)'); %xlim([0, xmax])
        figure(f(4));
        hs4 = plot(result.SA4(3,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of ASD');
        xlabel('ASD (degree)');ylabel('CDF (%)'); xlim([0, 120])
        figure(f(5));
        hs5 = plot(result.SA4(1,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of ZSD');
        xlabel('ZSD (degree)');ylabel('CDF (%)'); xlim([0, 50])
        figure(f(6));
        hs6 = plot(result.SA4(4,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of ASA');
        xlabel('ASA (degree)');ylabel('CDF (%)'); xlim([0, 120])
        figure(f(7));
        hs7 = plot(result.SA4(2,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(jj),...
            'LineWidth',1.5,...
            'LineStyle','none','Color',color(jj,:));grid on;hold on;
%         title('CDF of ZSA');
        xlabel('ZSA (degree)');ylabel('CDF (%)'); xlim([0, 50])
        hs = [hs,hs1,hs2,hs3,hs4,hs5,hs6,hs7];
    end
end
%%
TR = '38901';
freq = [6e9, 3e10, 6e10, 7e10];
lins = '- -.--: ';
for ii = 3
    for jj = 1:4
        frequency = freq(jj);
        if ii == 1
            name = '3D-UMa'; xmax = 2000; xm = -40;
        elseif ii == 2
            name = '3D-UMi'; xmax = 2000; xm = -40;
        elseif ii == 3
            name = '3D-InH'; xmax = 200;  xm = -30;
        end
        fname = ['.\results\TR',TR,'\full_config1_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
%         figure(f(((ii-1)*4+jj-1)*7+1));
        [~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
        figure(f(1));
        plot(result.CouplingLoss, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('Coupling loss (dB)');ylabel('CDF (%)'); 
%         if xm == -30
%             xlim([-100, xm])
%         else
%             xlim([floor(min(result.CouplingLoss)/10)*10, xm])
%         end
%         figure(f(((ii-1)*4+jj-1)*7+2));
        figure(f(2));
        plot(result.SIR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF (%)'); % xlim([floor(min(result.SIR)/10)*10, 30])
%         figure(f(((ii-1)*4+jj-1)*7+3));
        figure(f((3)));
        plot(result.DS, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of DS');
        xlabel('DS (ns)');ylabel('CDF (%)'); %xlim([0, xmax])
        figure(f(4));
        plot(result.SA4(3,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of ASD');
        xlabel('ASD (degree)');ylabel('CDF (%)'); xlim([0, 120])
        figure(f(5));
        plot(result.SA4(1,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of ZSD');
        xlabel('ZSD (degree)');ylabel('CDF (%)'); xlim([0, 50])
        figure(f(6));
        plot(result.SA4(4,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of ASA');
        xlabel('ASA (degree)');ylabel('CDF (%)'); xlim([0, 120])
        figure(f(7));
        plot(result.SA4(2,:), result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',size(jj),...
            'Marker',marker(2),...
            'LineWidth',2,...
            'LineStyle',lins(2*jj-1:2*jj),'Color',color(jj,:));grid on;hold on;
%         title('CDF of ZSA');
        xlabel('ZSA (degree)');ylabel('CDF (%)'); xlim([0, 50])
    end
end

%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    ylim([0 100]);
% m7 = plot(0,-10,'k-',...
%     'LineWidth',1.5); hold on;
% m8 = plot(0,-10,'k--s',...
%     'MarkerSize',11,...
%     'LineWidth',1.5); hold on;
% m9 = plot(0,-10,'k-.x',...
%     'MarkerSize',11,...
%     'LineWidth',1.5);
% 
% m1 = scatter(0,-10,'k',...
%     'Marker','o',...
%     'LineWidth',2); hold on;
% m2 = scatter(0,-10,'k',...
%     'Marker','^',...
%     'LineWidth',2); hold on;
% m3 = scatter(0,-10,'k',...
%     'Marker','s',...
%     'LineWidth',2);
% 
% m4 = plot(0,-10,'-',...
%     'LineWidth',3,...
%     'Color',[0 0.447058823529412 0.741176470588235]); hold on;
% m5 = plot(0,-10,'-',...
%     'LineWidth',3,...
%     'Color',[0.850980392156863 0.325490196078431 0.098039215686274]); hold on;
% m6 = plot(0,-10,'-',...
%     'LineWidth',3,...
%     'Color',[0.4 0.6 0.12]);
% 
% m10 = plot(0,-10,'-',...
%     'LineWidth',3,...
%     'Color',[0.494117647058824 0.184313725490196 0.556862745098039]);

legend([hr(a:7:end), hs(a:7:end)],'6 GHz, ref. [49]','30 GHz, ref. [49]','60 GHz, ref. [49]','70 GHz, ref. [49]','6 GHz, sim.','30 GHz, sim.','60 GHz, sim.','70 GHz, sim.','fontsize',14)
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
%     legend('UMa','UMa, mean of R1-165975','UMa, R1-165975','location','southeast');
%     saveas(f(a),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMa_',num2str(a),'.png']);
%     
%     figure(f(a+28));
%     legend('UMi','UMi, mean of R1-165975','UMi, R1-165975','location','southeast');
%     saveas(f(a+28),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
%     
%     figure(f(a+56));
%     legend('InH','InH, mean of R1-165975','InH, R1-165975','location','southeast');
%     saveas(f(a+56),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_InH_',num2str(a),'.png']);
end



%%
for a = 1:numel(picname)
    figure(f(a));
    ylim([0 100]);

m1 = plot(0,-13,'-x',...
    'LineWidth',2,'MarkerSize',11,...
    'Color',[0 0.447058823529412 0.741176470588235]); hold on;
m2 = plot(0,-15,'-.x',...
    'LineWidth',2,'MarkerSize',11,...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274]); hold on;
m3 = plot(0,-17,'--x',...
    'LineWidth',2,'MarkerSize',11,...
    'Color',[0.4 0.6 0.12]);

m4 = plot(0,-19,':x','MarkerSize',11,...
    'LineWidth',2,...
    'Color',[0.494117647058824 0.184313725490196 0.556862745098039]);

legend([m4,m3,m2,m1],'70 GHz','60 GHz','30 GHz','6 GHz','fontsize',14)
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
%     legend('UMa','UMa, mean of R1-165975','UMa, R1-165975','location','southeast');
%     saveas(f(a),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMa_',num2str(a),'.png']);
%     
%     figure(f(a+28));
%     legend('UMi','UMi, mean of R1-165975','UMi, R1-165975','location','southeast');
%     saveas(f(a+28),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
%     
%     figure(f(a+56));
%     legend('InH','InH, mean of R1-165975','InH, R1-165975','location','southeast');
%     saveas(f(a+56),['.\results\TR',TR,'\full_config1_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_InH_',num2str(a),'.png']);
end