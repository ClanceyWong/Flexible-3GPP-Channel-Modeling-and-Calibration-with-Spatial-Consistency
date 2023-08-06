function [UE_list, ue_pos_list] = Drop_UE(BSsector_list, antennaParams, ue_per_sec, R, min_d, beta, varargin)
    UE_list = []; ue_pos_list = [];
    ang.beta = beta;
    ang.gamma = 0;
    lowhigh = {'low','high'};
    id      = 0; 
    if ~isempty(varargin) && length(varargin) == 2
        num_din = length(BSsector_list)/3;
    else
        num_din = 1;
    end
    for uk = 1:ue_per_sec
        for m = 1:numel(BSsector_list)
            Ue      = elements.UE();
            id      = id + 1;
            Ue.ID   = id;
            fc      = BSsector_list(m).attached_BS.frequency;
            BW      = BSsector_list(m).BW;
            Ue.fcin = (BW/20)*(floor(20*rand)-(19/2)) + fc;
            if rand < 0.8
                Ue.bIndoor  = true;
                Ue.rand_din = 25*rand(1,num_din); % ,numel(BSsector_list)/3
%                 Ue.pos      = network_layout.drop_in_hexagon(BSsector_list(m).attached_BS.Position, R, min_d+Ue.rand_din(ceil(m/3)), BSsector_list(m).boresight, 360/numel(BSsector_list(m).attached_BS.sector));
                if num_din>1
                    Ue.pos      = network_layout.drop_in_hexagon2(BSsector_list(m).attached_BS.Position, R, min_d+Ue.rand_din(ceil(m/3)), BSsector_list(m).boresight);
                else
                    Ue.pos      = network_layout.drop_in_hexagon2(BSsector_list(m).attached_BS.Position, R, min_d+Ue.rand_din(1), BSsector_list(m).boresight);
                end
                N_fl        = randi([4, 8],1);
                Ue.n_fl     = randi([1, N_fl],1);
                Ue.O2IPL    = lowhigh{round(rand)+1};
            else
                Ue.rand_din = zeros(1,num_din);
%                 Ue.pos      = network_layout.drop_in_hexagon(BSsector_list(m).attached_BS.Position, R, (min_d), BSsector_list(m).boresight, 360/numel(BSsector_list(m).attached_BS.sector));
                Ue.pos      = network_layout.drop_in_hexagon2(BSsector_list(m).attached_BS.Position, R, min_d, BSsector_list(m).boresight);
                Ue.n_fl     = 1;
            end
            Ue.O2Isigma = randn;
            Ue.carPL = [9,5]; % for metallized car windows, ¦Ì = 20 can be used.
            Ue.rand_los = rand(1,numel(BSsector_list)/3);
            Ue.h_UT     = 3*(Ue.n_fl-1)+1.5;
            Ue.initPos  = [Ue.pos, Ue.h_UT].';
            ang.alpha   = rand*360-180;
%             panel       = antennas.antenna_panel(antennaParams.panel);
%             Ue.antenna  = antennas.antenna_array(antennaParams.array, ang, panel);
            Ue.antenna  = antennas.antenna_array(antennaParams, ang);
            Ue.antenna.attachedDevice = Ue;
            Ue.antenna.attachedType   = 'UE';
            ue_pos_list = [ue_pos_list; Ue.pos, Ue.h_UT]; %#ok<AGROW>
            UE_list     = [UE_list;Ue]; %#ok<AGROW>
        end
    end
    if ~isempty(varargin) && varargin{1}
        figure(1);
        for i = 1:numel(BSsector_list)
            indx = i:numel(BSsector_list):numel(UE_list);
            plot3(ue_pos_list(indx,1),ue_pos_list(indx,2),ue_pos_list(indx,3),'.','color',BSsector_list(i).lineColor);hold on;
        end
    end
    axis equal;view(0,90);
end

