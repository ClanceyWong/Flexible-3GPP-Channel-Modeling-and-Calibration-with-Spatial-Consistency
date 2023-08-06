function [azimuth,elevation,points] = fibonacci_sphere(samples)

    phi    = pi * (3. - sqrt(5.));          % golden angle in radians  0.382
    i      = 0:samples-1;
    y      = 1 - (i / (samples - 1)) .* 2;  % y goes from 1 to -1
    radius = sqrt(1 - y .* y);              % radius at y
    theta  = phi .* i;                      % golden angle increment
    x      = cos(theta) .* radius;
    z      = sin(theta) .* radius;
    [azimuth,elevation,~] = cart2sph(x,y,z);
    azimuth   = azimuth(:);
    elevation = elevation(:);
    points    = [x;y;z]';
end