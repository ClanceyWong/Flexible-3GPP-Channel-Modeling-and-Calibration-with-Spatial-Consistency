function pos = drop_in_office(bs_pos, min_d, xy_range)
    d = 0;
    while d < size(bs_pos,1)
        pos = [(xy_range(1,2)-xy_range(1,1))*rand+xy_range(1,1),(xy_range(2,2)-xy_range(2,1))*rand+xy_range(2,1)];
        pos_tmp = repmat(pos,size(bs_pos,1),1);
        diff = pos_tmp - bs_pos;
        dd = sqrt((diff(:,1)).^2 + (diff(:,2)).^2);
        d = sum(dd>min_d);
    end
end