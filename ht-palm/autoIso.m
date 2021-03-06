function isoVal = autoIso(vol)
%automatic calculation of iso-surface value for 
%volumetric data using otsu's method

vol= double(vol);
nBins = freedmanDiaconis(vol(:));
[n,xout] = hist(vol(:),nBins);

%from the multiOtsu m-file
[histo,pixval] = hist(vol(:),nBins);
P = histo/sum(histo);

%% Zeroth- and first-order cumulative moments
w = cumsum(P);
mu = cumsum((1:nBins).*P);

%% Maximal sigmaB^2 and Segmented image
sigma2B =...
    (mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));
[maxsig,k] = max(sigma2B);

isoVal = pixval(k+1);
