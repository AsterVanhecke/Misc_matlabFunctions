function mTrackErrorEst= getMTrackErrorEstimate(p,pcCell,mesh,fitquality,roiImg,roiExtDx,roiExtDy,roiImThresh,roiBox)

nTrim = 10;
nGradSmooth = 20;

if nargin == 0 || numel(pcCell) ==0
   mTrackErrorEst.imForceTot=0;
   mTrackErrorEst.imForceNorm=0;
   mTrackErrorEst.ftq=0;
   mTrackErrorEst.nPixInOut=0;
   mTrackErrorEst.nPixInOutNorm=0;
   mTrackErrorEst.length=0;
   mTrackErrorEst.width=0;
   mTrackErrorEst.wiggle=0;
   mTrackErrorEst.wiggleNorm=0;
else
   mTrackErrorEst.fitquality= fitquality;
   [mTrackErrorEst.imForceTot,mTrackErrorEst.imForceNorm] = getImageForce(p,pcCell,roiExtDx,roiExtDy,roiBox);
   [mTrackErrorEst.nPixInOut,mTrackErrorEst.nPixInOutNorm] = getPixInOut(p,pcCell,roiImThresh,roiBox);
   [mTrackErrorEst.length,mTrackErrorEst.width,mTrackErrorEst.wiggle,mTrackErrorEst.wiggleNorm] = getMeshInfo(mesh,nTrim ,nGradSmooth);
end

%-------------------------------------------------------------
function [imForceTot,imForceNorm] =getImageForce(p,pcCell,roiExtDx,roiExtDy,roiBox);
% Try to estimate quality of fit based on residual image forces

xCell=pcCell(:,1);
yCell=pcCell(:,2);
% Vector image forces (Fix,Fiy)
Fix =-p.imageforce * interp2a(1:roiBox(3)+1,1:roiBox(4)+1,roiExtDx,xCell,yCell,'linear',0);
Fiy = p.imageforce * interp2a(1:roiBox(3)+1,1:roiBox(4)+1,roiExtDy,xCell,yCell,'linear',0);

imForceTot =  sum(sqrt(Fix.^2+Fiy.^2))./p.imageforce;
imForceNorm = imForceTot/length(Fix)/2;

%-------------------------------------------------------------
function [nPixInOut,nPixInOutNorm] = getPixInOut(p,pcCell,roiImThresh,roiBox)
% as per Ulmann et al

nPixInOut=0;
nPixInOutNorm=0;
xCell=pcCell(:,1);
yCell=pcCell(:,2);

%convert cell outline to mask
fittedCell = poly2mask(xCell,yCell,roiBox(4)+1,roiBox(3)+1);
differenceIm = (roiImThresh~=fittedCell);
nPixInOut = sum(differenceIm(:));
nPixInOutNorm = nPixInOut/sum(fittedCell(:));

%-----------------------------------------------------------------
function [length,maxWidth,wiggle,wiggleNorm] = getMeshInfo(mesh,nTrim ,nGradSmooth)

centreLine = [mean([mesh(:,1),mesh(:,3)],2), mean([mesh(:,2),mesh(:,4)],2)];

dx = diff(centreLine(:,1));
dy = diff(centreLine(:,2));

%LENGTH
%length = sqrt(sum(dx.^2+dy.^2));
%BUG FIXED 131217!
length = sum(sqrt(dx.^2+dy.^2));

%WIDTH
width = sqrt( (mesh(:,1) - mesh(:,3)).^2 + (mesh(:,1) - mesh(:,3)).^2); 
maxWidth =max(width);

%WIGGLE
centreLineTrim = centreLine(nTrim+1:end-nTrim,:);
nPts = numel(centreLineTrim);
centreLineSmth = [smooth(centreLineTrim(:,1),nGradSmooth),smooth(centreLineTrim(:,2),nGradSmooth)];
dxSmth = diff(centreLineSmth(:,1));
dySmth = diff(centreLineSmth(:,2));
theta = atan2(dySmth,dxSmth)*180/pi;
dTheta = diff(theta);
wiggle = sum(abs(dTheta));
wiggleNorm = wiggle/nPts;



