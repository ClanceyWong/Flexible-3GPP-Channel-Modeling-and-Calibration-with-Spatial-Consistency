function [BS,BSsector_list] = Drop_BaseStation(scenario, antennaParams, sec_per_BS, frequency, version, varargin)
    %DROP_BASESTATION 此处显示有关此函数的摘要
    %   此处显示详细说明
    % drop BS 
    BS_pos_list = scenario.BSposition;
    if ~isempty(varargin) && varargin{1}
        figure(1); hold off;
        plot3(BS_pos_list(:,1),BS_pos_list(:,2),scenario.BS_height*ones(size(BS_pos_list,1),1),'ro','markersize',5,'linewidth',2);hold on;
        if strcmp(scenario.name,'3D-RMa') ||strcmp(scenario.name,'3D-UMa')||strcmp(scenario.name,'3D-UMi')
%             around = scenario.R*[1,0;0.5,0.866025403784439;-0.5,0.866025403784439;-1,0;-0.5,-0.866025403784439;0.5,-0.866025403784439;1,0];
            around2 = scenario.R/sqrt(3)*[sqrt(3)/2,-0.5;0,0;0,1;sqrt(3)/2,1.5;sqrt(3),1;sqrt(3),0;sqrt(3)/2,-0.5;sqrt(3)/2,-1.5;0,-2;-sqrt(3)/2,-1.5;-sqrt(3)/2,-0.5;-sqrt(3),0;-sqrt(3),1;-sqrt(3)/2,1.5;0,1;0,0;-sqrt(3)/2,-0.5];
            for i = 1:size(BS_pos_list,1)
                bspos = BS_pos_list(i,:) + around2;
                plot(bspos(:,1),bspos(:,2),'k');hold on;
            end
        elseif strcmp(scenario.name,'3D-InH')
            around = [-60, 25;60,25;60,-25;-60,-25;-60,25];
            plot(around(:,1),around(:,2),'k','linewidth',2); grid on;
        elseif strcmp(scenario.name,'InF')
            around = [scenario.x_range(1), scenario.y_range(1);scenario.x_range(2), scenario.y_range(1);scenario.x_range(2), scenario.y_range(2);scenario.x_range(1), scenario.y_range(2);scenario.x_range(1), scenario.y_range(1)];
            plot(around(:,1),around(:,2),'k','linewidth',2); grid on;
        end
        if length(varargin)==2 && strcmp(varargin{2},'SpatialConsistency')
            spatialconsist = true;
        elseif length(varargin)==2 && strcmp(varargin{2},'blockage')
            blockage = addition_components.blockage(frequency,'Landscape',scenario.name);
        end
    end
    BSsector_list = [];
    total         = 0;
    ang.beta  = 0; 
    ang.gamma = 0;
    for pos_idx = 1:size(BS_pos_list,1)
        BS_pos                = BS_pos_list(pos_idx,:);
        BS(pos_idx)           = elements.BaseStation(scenario); %#ok<*AGROW>
        BS(pos_idx).ID        = pos_idx;
        BS(pos_idx).initalPOS = BS_pos;
        BS(pos_idx).Position  = BS_pos;
        BS(pos_idx).frequency = frequency;
        if exist('blockage','var')
            BS(pos_idx).blockage  = blockage;
        end
        BS(pos_idx).getRandnGrid(version);
        if exist('spatialconsist','var') && spatialconsist
            BS(pos_idx).getSpatialConsistency();
        end
        for n = 1:sec_per_BS
            total               = total + 1;
            sector              = elements.BS_Sector(BS(pos_idx));
            sector.ID           = [pos_idx, n, total];
            BS(pos_idx).sector  = [BS(pos_idx).sector; total];
            sector.frequency    = frequency;
            sector.boresight    = (-90+(360/sec_per_BS)*(n-1));
            sector.BW           = scenario.BW;
            ang.alpha           = sector.boresight; 
%             panel               = antennas.antenna_panel(antennaParams.panel);
%             sector.antenna      = antennas.antenna_array(antennaParams.array, ang, panel);
            sector.antenna  = antennas.antenna_array(antennaParams, ang);
            sector.antenna.attachedDevice = sector;
            sector.antenna.attachedType   = 'BS';
            if strcmp(scenario.name,'InF')
                sector.boresight    = 0;
                sector.antenna.attachedType   = 'BSInF';
            end
            BSsector_list       = [BSsector_list; sector];
        end
    end
%     s1 = BS(1).corr_sos; len = length(BS(1).corr_sos);
%     for pos_idx = 2:size(BS_pos_list,1)
%         rn = 1+(round(rand(size(s1{1}.sos_phase,1),1))*2-1)*0.5;
%         if len == 3
%             BS(pos_idx).corr_sos{len}.sos_phase(:,1) = s1{len}.sos_phase(:,1).*rn(:,1);
%         end
%     end
%     BSsector_list = reshape(BSsector_list,3,[]);
end

