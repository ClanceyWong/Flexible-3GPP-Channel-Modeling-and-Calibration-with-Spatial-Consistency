function [ theta, phi, position ] = SpherePack( N )
    % PACK_SPHERE Creates equally distributed points on the unit sphere
    %   N:         The number of points to place.
    %   theta:     Vector of elevation angles in [rad] having values from -pi/2 to pi/2; size [ M x 1 ]
    %   phi:       Vector of azimuth angles in [rad] having values from -pi to pi; size [ M x 1 ]
    %   position:  Positions of the points on Cartesian coordinates; size [ 3 x M ]

    a       = 4*pi/N;           % Sphere surface
    d       = sqrt(a);          % Distance between points
    M_theta = round( pi/d );    % No. Latitudes
    d_theta = pi / M_theta;     % Resolution
    d_phi   = a / d_theta;
    A       = zeros(N+50,2);    % The output coordinates
    N_count = 1;                % Number of placed points
    
    for m = 1 : M_theta
        theta = pi * ( m-0.5 ) / M_theta;           % Current latitude
        M_phi = round( 2*pi*sin(theta)/d_phi );     % No. points on current latitude
        phi   = 2*pi*(0:M_phi-1) / M_phi;           % Longitude points
        N_new = numel( phi );
        
        A( N_count : N_count + N_new - 1 , 1 ) = theta;
        A( N_count : N_count + N_new - 1 , 2 ) = phi;
        N_count = N_count + N_new;
    end
    N_count = N_count - 1;

    % Format output
    if N_count < N && (N-N_count) <= 2
        A(N_count+1,:) = 0;
        A(N_count+2,:) = pi;
    elseif N_count < N && (N-N_count) > 2
        A(N_count+1,:)     = 0;
        A(N_count+2,:)     = pi;
        A((N_count+3):N,:) = [pi*rand((N-N_count),1);2*pi*rand((N-N_count),1)];
    elseif N_count > N
        [~,ind]        = sort(rand(N_count,1));
        A(1:N_count,:) = A(ind,:);
    end
    theta = A(1:N,1)-pi/2;
    phi   = angle( exp( 1j.*A(1:N,2) ) );

    if nargout > 2
        position = zeros( 3,N );
        position(1,:) = cos(theta).*cos(phi);
        position(2,:) = cos(theta).*sin(phi);
        position(3,:) = sin(theta);
    end
end
