function S = map( obj, xc, yc, zc )
    % MAP Generates a map at the given coordinates in x,y and z direction
    % This method generates a multi-dimensional array (of up to 3 dimensions) of spatially correlate random variables.  
    %
    % Input:
    %   xc      Vector of x coordinates in [m]
    %   yc      Vector of y coordinates in [m]
    %   zc      Vector of z coordinates in [m]
    %
    % Output:
    %   S       Array of spatially correlated random variables
    ndim_out = nargin-1;
    if ndim_out > obj.dimensions
        error( 'Number of requested dimensions is not supported.' );
    end

    if ~exist( 'yc','var' )
        yc = 0;
    end
    if ~exist( 'zc','var' )
        zc = 0;
    end

    nx = numel(xc);
    ny = numel(yc);
    nz = numel(zc);
    ox =  ones(nx,1,'uint8');
    oy =  ones(ny,1,'uint8');
    oz =  ones(nz,1,'uint8');

    x = reshape( single(xc) , 1, [] );
    x = x( oy,:,oz );
    x = x(:);
    y = reshape( single(yc) , [] , 1 );
    y = y( :,ox,oz );
    y = y(:);
    z = reshape( single(zc) , 1 , 1, []  );
    z = z( oy,ox,: );
    z = z(:);

    switch ndim_out
        case 1
            S = obj.value( x.' );
        case 2
            S = obj.value( [x,y].' );
        case 3
            S = obj.value( [x,y,z].' );
    end
    S = reshape( S, ny, nx, nz );

end
