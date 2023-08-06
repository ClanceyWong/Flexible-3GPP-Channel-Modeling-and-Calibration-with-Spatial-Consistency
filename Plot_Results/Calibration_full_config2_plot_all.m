clear
cd ..
TR = '38900';
date = '2020-10-16';
if ~exist(['.\results\TR',TR,'\full_config2_result_',date] ,'file')
    error("The result file doesn't exist.");
end

set(0,'defaultAxesFontName','Times')                                    % Default Font Type
set(0,'defaultTextFontName','Times')                                    % Default Font Type

picname = [{'CouplingLoss'},{ 'Geometry_SIR'},{'Largest_Singular_Value'},{'Smallest_Singular_Value'},{'Ratio_Singular_Value'}];
sheetName = ["UMa-6GHz","UMa-30GHz","UMa-60GHz","UMa-70GHz","UMi-6GHz","UMi-30GHz","UMi-60GHz","UMi-70GHz","InH-6GHz","InH-30GHz","InH-60GHz","InH-70GHz"];
f = [];
for fign = 1:numel(sheetName)*numel(picname)
    f = [f,figure];
end

%%
freq = [6e9, 3e10, 6e10, 7e10];
for ii = 1:3
    for jj = 1:4
        frequency = freq(jj);
        if ii == 1
            name = '3D-UMa'; xmax = 2000; xm = -40;xm2 = [-40,20];
        elseif ii == 2
            name = '3D-UMi'; xmax = 2000; xm = -40;xm2 = [-40,20];
        elseif ii == 3
            name = '3D-InH'; xmax = 200;  xm = -30;xm2 = [-50,10];
        end
        fname = ['.\results\TR',TR,'\full_config2_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
        figure(f(((ii-1)*4+jj-1)*5+1));
        plot(result.CouplingLoss, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('CouplingLoss (dB)');ylabel('CDF'); 
        if xm == -30
            xlim([-100, xm])
        else
            xlim([floor(min(result.CouplingLoss)/10)*10, xm])
        end
        figure(f(((ii-1)*4+jj-1)*5+2));
        plot(result.SIR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF'); xlim([floor(min(result.SIR)/10)*10, 30])
        figure(f(((ii-1)*4+jj-1)*5+3));
        plot(result.SV(2,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of largest singular value');
        xlabel('10log10(largest singular value)');ylabel('CDF'); xlim([-5,25]);
        figure(f(((ii-1)*4+jj-1)*5+4));
        plot(result.SV(1,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of smallest singular value');
        xlabel('10log10(smallest singular value)');ylabel('CDF'); xlim(xm2);
        figure(f(((ii-1)*4+jj-1)*5+5));
        plot(result.RSV, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ratio of singular value');
        xlabel('10log10(ratio singular value)');ylabel('CDF'); xlim([0,60]);
    end
end

for k = 1:12
    data3GPP = data.importfile(".\Docs\R1-165975_full_calibraton\Phase2Config2Calibration_v28_CMCC.xlsx", sheetName(k));
    for kk = 1:numel(picname)
        figure(f((k-1)*numel(picname)+kk));
%         figure(f(kk));
        plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',2); hold on;
        plot(data3GPP(:,((kk-1)*20+1):(kk*20-1)),0:100,'-','color',[0.7,0.7,0.7],'linewidth',2); hold on;
        plot(data3GPP(:,(kk*20)),0:100,'r-','linewidth',2); hold on;
    end
end
%%
for ii = 1:3
    for jj = 1:4
        frequency = freq(jj);
        if ii == 1
            name = '3D-UMa'; xmax = 2000; xm = -40;xm2 = [-40,20];
        elseif ii == 2
            name = '3D-UMi'; xmax = 2000; xm = -40;xm2 = [-40,20];
        elseif ii == 3
            name = '3D-InH'; xmax = 200;  xm = -30;xm2 = [-50,10];
        end
        fname = ['.\results\TR',TR,'\full_config2_result_',date,'\',name,'_',num2str(frequency/1e9),'GHz.mat'];
        load(fname,'result');
        figure(f(((ii-1)*4+jj-1)*5+1));
        plot(result.CouplingLoss, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of CouplingLoss');
        xlabel('CouplingLoss (dB)');ylabel('CDF'); 
        if xm == -30
            xlim([-100, xm])
        else
            xlim([floor(min(result.CouplingLoss)/10)*10, xm])
        end
        figure(f(((ii-1)*4+jj-1)*5+2));
        plot(result.SIR, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of Geometry SIR');
        xlabel('Geometry SIR (dB)');ylabel('CDF'); xlim([floor(min(result.SIR)/10)*10, 30])
        figure(f(((ii-1)*4+jj-1)*5+3));
        plot(result.SV(2,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of largest singular value');
        xlabel('10log10(largest singular value)');ylabel('CDF'); xlim([-5,25]);
        figure(f(((ii-1)*4+jj-1)*5+4));
        plot(result.SV(1,:), result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of smallest singular value');
        xlabel('10log10(smallest singular value)');ylabel('CDF'); xlim(xm2);
        figure(f(((ii-1)*4+jj-1)*5+5));
        plot(result.RSV, result.Pr,'b-', 'linewidth',2);grid on;hold on;
%         title('CDF of ratio of singular value');
        xlabel('10log10(ratio singular value)');ylabel('CDF'); xlim([0,60]);
    end
end
%%
% loc = {'southeast','northwest','southeast'};
for a = 1:numel(picname)*4
    figure(f(a));
    legend('UMa','UMa, mean of R1-165975','UMa, R1-165975','location','southeast');
    saveas(f(a),['.\results\TR',TR,'\full_config2_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMa_',num2str(a),'.png']);
    
    figure(f(a+20));
    legend('UMi','UMi, mean of R1-165975','UMi, R1-165975','location','southeast');
    saveas(f(a+20),['.\results\TR',TR,'\full_config2_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_UMi_',num2str(a),'.png']);
    
    figure(f(a+40));
    legend('InH','InH, mean of R1-165975','InH, R1-165975','location','southeast');
    saveas(f(a+40),['.\results\TR',TR,'\full_config2_result_',date,'\',cell2mat(picname(1+mod(a-1,numel(picname)))),'_InH_',num2str(a),'.png']);
end