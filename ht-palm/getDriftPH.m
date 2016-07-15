function [xDriftPix yDriftPix] = getDriftPH(phPre,phPost); 

phPre = double(phPre);
phPost= double(phPost);
[driftN, corAmplitudeN] = getImDrift(phPre,phPost);
xDriftPix = driftN(1);
yDriftPix = driftN(2);
%--------------------------------------------
function [drift, corAmplitude] = getImDrift(templateIm,featureIm, windowSize)
% spatial cross correlation function = SCCF
if ~exist('windowSize','var')
  windowSize= 10; %seems like a pretty reasonable number
end

%calculate the correlation function
% C is the image of the correlation function. 
%zeroCoord is the [i,j] coordinate corresponding to zero displacement in
% the correlation function 
[C,zeroCoord] = corrfunc(templateIm,featureIm);


% extract the drift by fitting a gaussian to the SCCF
% careful with (i,j) vs (x,y)!!

[corrMaxPos corAmplitude] =   getPosGauss2d(C,zeroCoord,windowSize);
%[corrMaxPos corAmplitude] =   getPeakPosCentroid(C,zeroCoord,windowSize); %this does not seem sufficiently accurate!
drift = zeroCoord  - corrMaxPos;

%%OPTIONAL PLOTTING OF THE CORRELATION IMAGES
%subplot(3,1,1)
%imagesc(templateIm);
%subplot(3,1,2)
%imagesc(featureIm);
%subplot(3,1,3)
%imagesc(C);

%%-------------------------------
%-------------------------------------------------------------
%function [G,zeroCoord] = getImDrift(template,feature)
function [G,zeroCoord] = corrfunc(template,feature)
% July 9, 2003
% David Kolin
% 18/03/11
% Seamus Holden
% 1)Minor modification 2011 SH to output zeroCoordinate in ICS image
% 2) 110318 Now a very heavily modified version of original function. Does cross not auto correlation
% NB:template and feature should be same size
% zeroCoord is (x,y) coordinates, not (i,j)!
template=double(template);
feature=double(feature); 

% Calculates 2D correlation functions for each z slice of a given 3D matrix (imgser)
% Output as a 3D matrix with same dimensions as imgser
G = zeros(size(template)); % Preallocates matrix

% Calculates corr func
%autocorrelation:
%%G = (fftshift(real(ifft2(fft2(template).*conj(fft2(template)))))) ...
%%        /(mean(template(:))^2*size(template,1)*size(template,2) ) ...
%%      -1;

%cross correlation
G = (fftshift(real(ifft2(fft2(template).*conj(fft2(feature))))))/...
		( (mean(template(:))*mean(feature(:))) * size(template,1)*size(template,2) ) ...
		- 1;
% SH mod
% make sure that you know where the zero coordinate is from the DFT
% so that we can calculate absolute drift
imsize = size(template);
zeroCoordX = (floor(imsize(2)/2)+1);
zeroCoordY = (floor(imsize(1)/2)+1);
zeroCoord = [ zeroCoordX,zeroCoordY];

%----------------------------------------------------------------------------------------------
