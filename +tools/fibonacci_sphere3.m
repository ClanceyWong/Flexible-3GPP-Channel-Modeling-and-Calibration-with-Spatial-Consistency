function [azimuth,elevation,points] = fibonacci_sphere3(samples)
 
    phi    = 2*pi * ((sqrt(5.)+1)/2 - 1);           % golden angle in radians 0.618
    i      = 1:samples;
    z      = (2*i-1) / samples - 1;       % y goes from 1 to -1
    radius = sqrt(1 - z .* z);              % radius at y
    theta  = phi .* i;                      % golden angle increment
    x      = cos(theta) .* radius;
    y      = sin(theta) .* radius;
    [azimuth,elevation,~] = cart2sph(x,y,z);
    azimuth   = azimuth(:);
    elevation = elevation(:);
    points    = [x;y;z]';
end