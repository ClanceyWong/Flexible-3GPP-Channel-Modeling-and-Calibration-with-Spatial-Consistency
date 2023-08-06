function [azimuth,elevation,points] = uniform_sphere(samples)
    
    theta     = 2*pi*rand(1,samples);
    phi       = acos(2*rand(1,samples)-1);
    x         = sin(phi).*cos(theta);
    y         = sin(phi).*sin(theta);
    z         = cos(phi);
    azimuth   = theta(:);
    elevation = phi(:);
    points    = [x;y;z]';
end