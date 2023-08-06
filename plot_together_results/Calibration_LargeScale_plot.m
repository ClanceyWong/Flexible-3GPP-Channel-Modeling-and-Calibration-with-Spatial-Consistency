clear
cd ..
TR   = '38900';
date = '2020-08-18';
if ~exist(['.\results\TR',TR,'\large_scale_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type
set(0,'defaultAxesFontSize',13)                                    % Default Font Type
set(0,'defaultTextFontSize',13)                                    % Default Font Type
f1 = figure('InvertHardcopy','off','Color',[1 1 1],'Renderer','painters',...
    'OuterPosition',[897 343 576 514]); 
f2 = figure('InvertHardcopy','off','Color',[1 1 1],'Renderer','painters',...
    'OuterPosition',[897 343 576 514]); 
f3 = figure('InvertHardcopy','off','Color',[1 1 1],'Renderer','painters',...
    'OuterPosition',[897 343 576 514]);
f = [f1,f2,f3];
picname = [{'CouplingLoss'},{ 'Geometry_SINR'},{ 'Geometry_SIR'}];
indx = [2,22;33,53;64,84];
hr = [];
 marker = 'so^d';
 msize = [12,12,12,12];
 color = [0 0.447058823529412 0.741176470588235;...
     0.850980392156863 0.325490196078431 0.098039215686274;...
     0.4 0.6 0.12;...
     0.494117647058824 0.184313725490196 0.556862745098039];

sheetName = ["indoor-6GHz","indoor-30GHz","indoor-70GHz"]; % "UMi-6GHz","UMi-30GHz","UMi-70GHz","UMa-6GHz","UMa-30GHz","UMa-70GHz",
for k = 1:length(sheetName)
    data3GPP = data.importfile_p1_all(".\Docs\R1-165974_large_scale_calibration\Phase1Calibration_v42_CMCC.xlsx", sheetName(k), [29, 129]);
    for kk = 1:numel(picname)
%         figure(f((ceil(k/3)-1)*3+kk));
        figure(f(kk))
        plot(str2double(data3GPP(:,indx(kk,1):indx(kk,2))),0:100,'color',[.85,.85,.85],'linewidth',3); hold on;
    end
end
for k = 1:length(sheetName)
    data3GPP = data.importfile_p1(".\Docs\R1-165974_large_scale_calibration\Phase1Calibration_v42_CMCC.xlsx", sheetName(k), [29, 129]);
    ind = floor((k-1)/3)+1;
    inc = mod(k-1,3)+1;
    for kk = 1:numel(picname)
%         figure(f((ceil(k/3)-1)*3+kk));
        figure(f(kk));
        hrr = plot(data3GPP(:,kk),0:100,'MarkerIndices',[5 15 25 35 45 55 65 75 85 95]+1,'MarkerSize',msize(ind),... %    'Marker',marker(ind),...
            'LineWidth',2,...
            'LineStyle','-','Color',color(inc,:)); hold on;
        hr = [hr,hrr];
    end
end
%%
hs900 = [];
frequency = [6e9, 3e10,7e10];
for jj = 1:3
    for ii = 3
        if ii == 2
            name = '3D-UMa';
        elseif ii == 1
            name = '3D-UMi'; 
        elseif ii == 3
            name = '3D-InH'; 
        end
        fname = ['.\results\TR38900\large_scale_result_',date,'\',name,'_',num2str(frequency(jj)/1e9),'GHz.mat'];
        load(fname,'result');
        [~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
%         figure(f((ii-1)*3+1));
        figure(f(1));
        hs1 = plot(result.CouplingLoss, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(ii),...
            'Marker',marker(ii),...
            'LineWidth',1.5,...
            'LineStyle','-.','Color',color(jj,:));grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('Coupling loss (dB)');ylabel('CDF (%)'); %xlim([floor(min(result.CouplingLoss)/10)*10, -30])
        ylim([0,100]);
%         figure(f((ii-1)*3+2));
        figure(f(2));
        hs2 = plot(result.SINR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(ii),...
            'Marker',marker(ii),...
            'LineWidth',1.5,...
            'LineStyle','--','Color',color(jj,:));grid on;hold on;
%         title('CDF of Geometry SINR');
        xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); %xlim([floor(min(result.SINR)/10)*10, 30])
        ylim([0,100]);
%         figure(f((ii-1)*3+3));
        figure(f(3));
        hs3 = plot(result.SIR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(ii),...
            'Marker',marker(ii),...
            'LineWidth',1.5,...
            'LineStyle','-.','Color',color(jj,:));grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF (%)'); %xlim([floor(min(result.SIR)/10)*10, 30])
        ylim([0,100]);
        hs900 = [hs900,hs1,hs2,hs3];
    end
end
%%
% marker = 'xso^d';
%  msize = [11,11,11,7];
% hs901 = [];
% frequency = [6e9, 3e10,7e10];
% for jj = 1:3
%     for ii = 1:1
%         if ii == 2
%             name = '3D-UMa';
%         elseif ii == 1
%             name = '3D-UMi'; 
%         elseif ii == 3
%             name = '3D-InH'; 
%         end
%         fname = ['.\results\TR38901\large_scale_result_',date,'\',name,'_',num2str(frequency(jj)/1e9),'GHz.mat'];
%         load(fname,'result');
%         [~,start] = min(abs(result.Pr-5)); step = numel(result.Pr)/10;
% %         figure(f((ii-1)*3+1));
%         figure(f(1));
%         hs1 = plot(result.CouplingLoss, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(ii),...
%             'Marker',marker(ii),...
%             'LineWidth',1.5,...
%             'LineStyle','-','Color',color(jj,:));grid on;hold on;
% %         title('CDF of CouplingLoss');
%         xlabel('Coupling loss (dB)');ylabel('CDF (%)'); %xlim([floor(min(result.CouplingLoss)/10)*10, -30])
%         ylim([0,100]);xlim([-220,-60])
% %         figure(f((ii-1)*3+2));
%         figure(f(2));
%         hs2 = plot(result.SINR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(ii),...
%             'Marker',marker(ii),...
%             'LineWidth',1.5,...
%             'LineStyle','-','Color',color(jj,:));grid on;hold on;
% %         title('CDF of Geometry SINR');
%         xlabel('Geometry SINR (dB)');ylabel('CDF (%)'); %xlim([floor(min(result.SINR)/10)*10, 30])
%         ylim([0,100]);xlim([-100,40])
% %         figure(f((ii-1)*3+3));
%         figure(f(3));
%         hs3 = plot(result.SIR, result.Pr,'MarkerIndices',(start:step:numel(result.Pr)),'MarkerSize',msize(ii),...
%             'Marker',marker(ii),...
%             'LineWidth',1.5,...
%             'LineStyle','-','Color',color(jj,:));grid on;hold on;
% %         title('CDF of Geometry SIR');
%         xlabel('Geometry SIR (dB)');ylabel('CDF (%)'); %xlim([floor(min(result.SIR)/10)*10, 30])
%         ylim([0,100]);
%         hs901 = [hs901,hs1,hs2,hs3];
%     end
% end

%%
loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
%     m7 = plot(0,-10,'k-',...
%     'LineWidth',1.5); hold on;
% m8 = plot(0,-10,'k--s',...
%     'MarkerSize',12,...
%     'LineWidth',1.5); hold on;
% m9 = plot(0,-10,'k-.x',...
%     'MarkerSize',12,...
%     'LineWidth',1.5);

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
% m4 = plot(0,-10,...
%     'LineWidth',3,...
%     'Color',[0 0.447058823529412 0.741176470588235]); hold on;
% m5 = plot(0,-10,...
%     'LineWidth',3,...
%     'Color',[0.850980392156863 0.325490196078431 0.098039215686274]); hold on;
% m6 = plot(0,-10,...
%     'LineWidth',3,...
%     'Color',[0.4 0.6 0.12]);
% m1,m2,m3,,'UMa','UMi','InH'
legend([hr(a:3:end), hs900(a:3:end)],'6 GHz, ref. [48]','30 GHz, ref. [48]','70 GHz, ref. [48]','6 GHz, sim.','30 GHz, sim.','70 GHz, sim.','fontsize',14)
grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
%     saveas(f(a+6),['.\results\TR',TR,'\large_scale_result_',date,'\',cell2mat(picname(a)),'_InH.png']);
end

%% only 38901

loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)
    figure(f(a));
    m7 = plot(0,-10,'k-',...
    'LineWidth',1.5); hold on;
m8 = plot(0,-10,'k--s',...
    'MarkerSize',12,...
    'LineWidth',1.5); hold on;
m9 = plot(0,-10,'k-.x',...
    'MarkerSize',12,...
    'LineWidth',1.5);

m1 = plot(0,-1,'k',...
    'Marker','s',...
    'LineWidth',2); hold on;
m2 = plot(0,-3,'k',...
    'Marker','x',...
    'LineWidth',2); hold on;
m3 = plot(0,-5,'k',...
    'Marker','o',...
    'LineWidth',2);

m4 = plot(0,-7,...
    'LineWidth',2,...
    'Color',[0 0.447058823529412 0.741176470588235]); hold on;
m5 = plot(0,-9,'-.',...
    'LineWidth',2,...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274]); hold on;
m6 = plot(0,-11,'--',...
    'LineWidth',2,...
    'Color',[0.4 0.6 0.12]);
% ,'UMa','UMi','InH'
legend([m1,m2,m3,m6,m5,m4], 'UMa','UMi','InH','6 GHz','30 GHz','60 GHz','fontsize',14)
grid on; set(gca,'GridLineStyle','--','GridColor',[0.75,0.75,0.75],'GridAlpha',1);
%     saveas(f(a+6),['.\results\TR',TR,'\large_scale_result_',date,'\',cell2mat(picname(a)),'_InH.png']);
end