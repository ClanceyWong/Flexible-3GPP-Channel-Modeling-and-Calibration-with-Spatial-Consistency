function [ fr, Rhs ] = iteration_search( R, D, al )
    % Detect a single multipath component using the SAGE algorithm
    %
    %   R           The ACF
    %   D           The normalized distances
    %   al          The desired amplitude
    %   fr          The sample frequency
    %   Rhs         The reconstructed ACF 

    % Global search
    fr = -15 : 0.1 : 15;                % The test-frequency range
    fr = single( fr );

    % Calculate the "cost-function" for each test-frequency
    x = al * cos( fr.' * D );
    x = R( ones(1,numel(fr)) ,:) - x;
    x = sum( x.^2, 2 );

    % There migt be 2 ideintcal values due to symmetry reasons, pick one randomly.
    [x,ind] = sort(x);
    if abs(x(2)-x(1)) < 1e-6
        ind = ind( randi(2) );
    else
        ind = ind(1);
    end
    x = x(1);

    % Start point an search direction
    fr = fr(ind);
    stp = single( 0.037 );

    while abs( stp ) > 1e-6
        % Update the frequency
        frn = fr + stp;

        % This implements the sum as a vector product
        xn = sum( ( R - al * cos( frn * D ) ).^2 );

        if xn >= x
            stp = -0.21 * stp;
            fr  = fr + stp;
        else
            x  = xn;
            fr = frn;
        end
    end

    % Update with the last value from the loop
    Rhs = al * cos( frn * D );
end