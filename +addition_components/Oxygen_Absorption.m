classdef Oxygen_Absorption
    % OXYGEN_ABSORPTION 
    
    % The final frequency-domain channel response is obtained by 
    % the summation of frequency-domain channel responses of all clusters.
    % Time-domain channel response is obtained by the reverse transform 
    % from the obtained frequency-domain channel response.
    
    
    properties
        % {tau_delta} is 0 in the LOS case and min(tau_n') otherwise, 
        % where min(tau_n') is the minimum delay in Step 5.
        tau_delta    
        alpha_fc
        LBW_enable = false
    end
    properties(Constant)
        c              = 3e8  % m/s
        freq           = [0, 52, 53,  54, 55,  56,  57,   58,   59, 60,   61,   62,   63,  64,  65,  66, 67, 68, 100]
        Loss_freq_dBkm = [0,  0,  1, 2.2,  4, 6.6, 9.7, 12.6, 14.6, 15, 14.6, 14.3, 10.5, 6.8, 3.9, 1.9,  1,  0,   0]
    end
    
    methods
        function obj = Oxygen_Absorption(fc, tau_delta, delta_f)
            % OXYGEN_ABSORPTION 
            obj.alpha_fc = obj.interp(obj.freq,[],obj.Loss_freq_dBkm,fc);
            obj.tau_delta = tau_delta;
            if exist('delta_f','var') % large BandWidth enable
%               For large channel bandwidth, first transform the time-domain channel response of each cluster 
%               (all rays within one cluster share common oxygen absorption loss for simplicity) into 
%               frequency-domain channel response, and apply the oxygen absorption loss to the cluster's 
%               frequency-domain channel response for frequency fc + ¦¤f within the considered bandwidth.
                f = fc + delta_f; % ¦¤f is in [-B/2, B/2], where B is the bandwidth.
                obj.alpha_fc = obj.interp(obj.freq,[],obj.Loss_freq_dBkm,f);
                obj.LBW_enable = true;
            end
        end
        
        function [lossdB, loss] = loss(obj,d3D,tau_n)
            % d3D [1x1]
            alpha_f = obj.alpha_fc(:);
            alpha_f = repmat(alpha_f,1,size(tau_n,2));
            if obj.LBW_enable % tau_n [1 x N]
                tau_n = repmat(tau_n(:).',size(alpha_f,1),1);
                lossdB = (alpha_f/1000)*(d3D + obj.c*(tau_n + obj.tau_delta));
            else % tau_n [size(alpha_f,1) x N]
                lossdB = (alpha_f/1000)*(d3D + obj.c*(tau_n + obj.tau_delta));
            end
            loss = 10.^(-lossdB/10);
        end
    end
    methods(Static)
        function zi = interp( x, y, z, xc, yc )
            %INTERP 2D linear interpolation optimized for performace
            % Description:
            %   This function implements a 2D linear interpolation which is highly optimized for fast execution. All calculations are done in single-precision floating point (30% faster than double
            %   precision, but less accurate), and multiple data sets can be interpolated simultaneously.  
            %   One-dimensional linear interpolation can be done by using:  zi = interp( x, 0, z, xc )
            if isempty( y )
                y = 0;
            end
            if isa(z,'double')
                z = single( z );
                z_is_double = true;
            else
                z_is_double = false;
            end
            x = single( x(:).' );
            y = single( y(:).' );
            nx = numel(x);
            ny = numel(y);
            ne = size( z,3 );
            if size( z,1 ) ~= ny || size( z,2 ) ~= nx
                error('Size of z does not match');
            end
            % Option for 1D linear interpolation
            if ny == 1
                y  = 0;
                yc = 0;
                z  = z(:);
            end
            nxc = numel(xc);
            nyc = numel(yc);
            oxc = ones(1,nxc,'uint8');
            oyc = ones(1,nyc,'uint8');
            xi = reshape( single(xc) , 1, [] );
            xi = xi( oyc,: );
            xi = xi(:).';
            yi = reshape( single(yc) , [] , 1 );
            yi = yi( :,oxc );
            yi = yi(:).';
            ni = numel(xi);
            ii = uint32( 1:ni );
            % Determine the nearest location of xi in x and the difference to
            % the next point
            [tmp,b] = sort( xi );
            [~,a]   = sort( [x,tmp] );
            ui      = uint32( 1:(nx + ni) );
            ui(a)   = ui;
            ui      = ui(nx+1:end) - ii;
            ui(b)   = ui;
            ui( ui==nx ) = nx-1;
            ui( ui==0 ) = 1;
            uin     = ui+1;
            u       = (xi-x(ui))./( x(uin)-x(ui) );
            u(isnan(u)) = 0;
            u       = u';
            % Determine the nearest location of yi in y and the difference to
            % the next point
            if ny > 1
                [tmp,b] = sort( yi );
                [~,a]   = sort( [y,tmp] );
                vi      = uint32( 1:(ny + ni) );
                vi(a)   = vi;
                vi      = vi(ny+1:end) - ii;
                vi(b)   = vi;
                vi( vi==ny ) = ny-1;
                vi( vi==0 ) = 1;
                vin     = vi+1;
                v       = (yi-y(vi))./( y(vin)-y(vi) );
                v(isnan(v)) = 0;
                v       = v';
            else
                vi  = uint32( 1 );
                vin = uint32( 1 );
                v   = zeros( ni,1,'single' );
            end
            % Calculate the scaling coefficients
            c1 = (1-v).*(1-u);
            c2 = (1-v).*u;
            c3 = v.*(1-u);
            c4 = v.*u;
            % Determine the indices of the elements
            pa = vi  + ( ui  -1 )*ny;
            pb = vi  + ( uin -1 )*ny;
            pc = vin + ( ui  -1 )*ny;
            pd = vin + ( uin -1 )*ny;
            pX = [pa,pb,pc,pd].';
            pY = uint32( (0:ne-1)*nx*ny );
            tr = true( ni,1 );
            fl = false( ni,1 );
            i1 = [tr;fl;fl;fl];
            i2 = [fl;tr;fl;fl];
            i3 = [fl;fl;tr;fl];
            i4 = [fl;fl;fl;tr];
            % Interpolate
            zi = zeros( ni, ne, 'single' );
            for n = 1 : ne
                ndx = pY(n) + pX;
                a = z( ndx );
                zi(:,n) = c1.*a(i1) + c2.*a(i2) + c3.*a(i3) + c4.*a(i4);
            end
            zi = reshape(zi,nyc,nxc,ne);
            if z_is_double
                zi = double( zi );
            end
        end
    end
end

