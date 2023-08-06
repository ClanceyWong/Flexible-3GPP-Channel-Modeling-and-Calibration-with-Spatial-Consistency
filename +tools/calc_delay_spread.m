function [ ds, mean_delay ] = calc_delay_spread(taus, pow )
% CALC_DELAY_SPREAD Calculates the delay spread in [s]
    N = size( taus,1 );
    if size( pow,1 ) < N
        pow = pow( ones(1,N), : );
    end
    % Normalize powers
    pt = sum( pow,2 );
    pow = pow./pt( :,ones(1,size(pow,2))  );

    mean_delay = sum( pow.*taus,2 ); 
    tmp = taus - mean_delay(:,ones( 1,size(taus,2) ) );
    ds = sqrt( sum(pow.*(tmp.^2),2) - sum( pow.*tmp,2).^2 );
end