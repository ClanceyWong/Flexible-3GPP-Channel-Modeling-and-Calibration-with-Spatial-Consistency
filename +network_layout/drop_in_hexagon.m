function pos = drop_in_hexagon(center, R, min_d, boresight, angle_range)
    d = 0; a = inf;
    while d <= min_d || a > (angle_range/2)
        v1 = R*exp(1j*pi/3)*rand;
        v2 = R*exp(-1j*pi/3)*rand;
        tmp = randsrc(1,1,[0, 1, 2]);
        y = (v1+v2)*exp(1j*2*pi/3*tmp);
        d = abs(y);
        deg = angle(y)*180/pi;
        a = min(abs(deg - boresight),abs(deg + 360 - boresight));
    end
    pos = [real(y), imag(y)] + center;
end