function xcorr = calc_xcorr(X,Y)
    if isreal(X) && isreal(Y)
        xcorr = (mean(X.*Y)-mean(X)*mean(Y))/sqrt(mean(X.^2)-mean(X)^2)/sqrt(mean(Y.^2)-mean(Y)^2);
    else
        xcorr = abs((mean(X.*conj(Y))-mean(X)*conj(mean(Y)))/sqrt(mean(X.*conj(X))-mean(X)*conj(mean(X)))/sqrt(mean(Y.*conj(Y))-mean(Y)*conj(mean(Y))));
    end
end

