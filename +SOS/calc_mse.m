function [ mse, mse_core, mse_all, Ro, Ri ] = calc_mse( sos, T )
% Description:
%   This method calculates the mean-square-error (in dB) of the approximation of the given ACF.
% Input:
%   T           Number of test directions, Default: 10
%
% Output:
%   mse         The MSE for the given number of test directions
%   mse_core    The 2D MSE for the range 0 to Dmax (covered by the desired ACF)
%   mse_all     The 2D MSE for the range 0 to 2*Dmax
%   Ro          Approximated ACF from sos.acf_approx
%   Ri          Interpolated ACF (from sos.acf_2d in case of 2 dimensions or more)

if ~exist( 'T','var' ) || isempty( T )
    T = 10;
end
if sos.dimensions == 1
    T = 1;
end

N = sos.no_coefficients;
S = numel( sos.dist(1):sos.dist(2):sos.dist(3) );
fr = sos.sos_freq;

if sos.dimensions == 1
    Ri = sos.acf;
    N = numel(Ri);
    Ri = [ zeros(1,N-1), Ri(end:-1:2), Ri, zeros(1,N-1) ];
    Ro = sos.acf_approx;
    
    mse_all  = -10*log10( mean( (Ro(:)-Ri(:)).^2 ) );
    ii       = 2*N-1 : 3*N-2;
    mse_core = -10*log10( mean( (Ro(ii)-Ri(ii)).^2 ) );
    mse = mse_core;
    
elseif sos.dimensions == 2
    
    % Generate test directions
    u = 2 * pi / T;
    phi_t = 0 : u : 2*pi;
    phi_t = single( angle(exp(1j*phi_t(1:T))) );
    
    % Test frequencies
    ft = fr(:,1) * cos( phi_t ) + fr(:,2) * sin( phi_t );
    
else
    
    % Generate test directions
    [ theta_t, phi_t ] = tools.SpherePack( T );
    theta_t = single( theta_t.' );
    phi_t = single( phi_t.' );
    T = numel( phi_t );
    
    % Test frequencies
    ft = fr(:,1) * ( cos( phi_t ) .* cos( theta_t ) ) +...
        fr(:,2) * ( sin( phi_t ) .* cos( theta_t ) ) +...
        fr(:,3) * sin( theta_t );
end

% Reconstructed ACF
Rh = ones( T, S, N, 'single' );
D = 2*pi*( sos.dist(1):sos.dist(2):sos.dist(3) );
oT = ones(1,T);

if sos.dimensions > 1
    
    for n = 1 : N
        Rh(:,:,n) = cos( ft(n,:)' * D );
    end
    dbg = sum( Rh,3 )/N;
    
    flag = 0;
    if flag 
        plot(sos.dist, dbg')
        hold on
        plot(( sos.dist(1):sos.dist(2):sos.dist(3) ),sos.acf,'-k','Linewidth',3)
        hold off
    end
    
    mse = sum( (sos.acf(oT,:) - dbg ).^2 ,2 ).' / S;
    mse = -10*log10( sum( mse ) / T );
    
    if nargout > 1
        
        Ri = sos.acf_2d;
        Ro = sos.acf_approx;
        
        N = numel(sos.acf);
        on = ones(1,4*N-3,'uint8');
        p = (1:4*N-3)-2*N+1;
        p = p.^2;
        ii = sqrt(p(on,:) + p(on,:)') <= N;
        
        if sos.dimensions == 2
            mse_all = -10*log10( mean( (Ro(:)-Ri(:)).^2 ) );
            mse_core = -10*log10( mean( (Ro(ii)-Ri(ii)).^2 ) );
        else
            for n = 1:3
                R3 = Ro(:,:,n);
                mse_all(n) = mean( (R3(:)-Ri(:)).^2 );
                mse_core(n) = mean( (R3(ii)-Ri(ii)).^2 );
            end
            mse_all = -10*log10( mean( mse_all ));
            mse_core = -10*log10( mean( mse_core ));
        end
    end
    
end

end

