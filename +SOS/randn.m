function [ s, obj ] = randn( dist_decorr, ca, cb, acf_type )
% Generates spatially correlated normal distributed random numbers
%
% Input:
%   dist_decorr:  Vector of decorrelation distances  [1 x M] or [ M x 1 ]
%   ca:           Coordinates for the first mobile device in [m] given as [3 x N] matrix. The rows correspond to the x,y and z coordinate.
%   cb:           Coordinates for the corresponding second mobile device in [m] given as [3 x N] matrix. This variable must either be empty or have the same size as "ca".
%   acf_type:     String describing the shape of the autocorrelation function and the number of sinusoids, Default: 'Comb300'
%
% Output:
%   s:            Random spatially correlated numbers [ M x N ]

s = zeros( numel( dist_decorr ), size(ca,2) );
if ~exist( 'cb','var' ) || isempty( cb )
   cb = []; 
end
if ~exist( 'acf_type','var' ) || isempty( acf_type )
   acf_type = 'Comb300'; 
end

obj = SOS.sos([]);
for n = 1 : numel( dist_decorr )
    if dist_decorr(n) == 0
        dist_decorr(n) = 0.1;
    end
    obj(1,n) = SOS.sos( acf_type, 'Normal', dist_decorr(n) );
    s(n,:) = obj(1,n).value( ca, cb );
end

end
