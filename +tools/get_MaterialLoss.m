function Mloss = get_MaterialLoss(fc)
    % fc in GHz
    %        glass   IRRglass   concrete  wood
    Mloss = [2+0.2*fc;23+0.3*fc;5+4*fc;4.85+0.12*fc];
end