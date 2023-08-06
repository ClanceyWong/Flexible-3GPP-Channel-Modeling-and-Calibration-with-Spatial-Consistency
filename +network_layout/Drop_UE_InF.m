function [UE_list, ue_pos_list] = Drop_UE_InF(BSsector_list, antennaParams, ue_per_sec, min_d, xy_range, beta, varargin)
    UE_list = []; ue_pos_list = [];
    ang.beta = beta;
    ang.gamma = 0;
    total = 0;
    bs = [BSsector_list(1:3:end).attached_BS];
    bspos = reshape([bs(:).Position]',2,[])';
    for uk = 1:ue_per_sec
        for m = 1:numel(BSsector_list)
            Ue          = elements.UE();
            total       = total + 1;
            Ue.ID       = total;
            fc          = BSsector_list(m).attached_BS.frequency;
            BW          = BSsector_list(m).BW;
            Ue.fcin     = (BW/20)*(floor(20*rand)-(19/2)) + fc;
            Ue.pos      = network_layout.drop_in_office(bspos, (min_d), xy_range);
            Ue.rand_los = rand(1,numel(BSsector_list));
            Ue.h_UT     = 1.5;
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