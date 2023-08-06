function pos = drop_in_hexagon2(center, R, min_d, boresight)
    d = 0;
    R = R/sqrt(3);
    center2 = center + R*[cosd(boresight),sind(boresight)];
    rota = [cosd(boresight),-sind(boresight);sind(boresight),cosd(boresight)];
    while d <= min_d
        v1 = R*exp(1j*pi/3)*rand;
        v2 = R*exp(-1j*pi/3)*rand;
        tmp = randsrc(1,1,[0, 1, 2]);
        y = (v1+v2)*exp(1j*2*pi/3*tmp);
        yt = (rota*[real(y); imag(y)])';
        pos = yt + center2;
        d = sqrt(sum((pos-center).^2));
    end
end

