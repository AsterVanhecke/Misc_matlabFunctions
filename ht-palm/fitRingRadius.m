function [c0,a,b] = fitRingRadius(x_th,y_rad,z,t,meshRot,varargin) 

radLim = [-500,500];
sigmaBlur = 15;
satVal = 0;
pixSize = 10;
gammaVal =1;
ampScaleFactor = 100;

ii = 1;
while ii <= numel(varargin)
  if strcmp(varargin{ii},'Lim')
    radLim = varargin{ii+1};
    ii = ii+2;
  elseif strcmp(varargin{ii},'Gamma')
    gammaVal = varargin{ii+1};
    ii = ii+2;
  elseif strcmp(varargin{ii},'SatVal')
    satVal= varargin{ii+1};
    ii = ii+2;
  elseif strcmp(varargin{ii},'BlurSigma')
    sigmaBlur= varargin{ii+1};
    ii = ii+2;
  elseif strcmp(varargin{ii},'PixSize')
    pixSize= varargin{ii+1};
    ii = ii+2;
  end
end

nCell = numel(t);
for ii = 1:nCell
%for ii = [1,2,3,19,21,26] %spot,ring, ring
  ii
  circleData = [y_rad{ii},z{ii}]';
  [c0(ii,:), a(ii), b(ii)] = fitEllipseOnAxisRobust2(circleData);

  %[c0, r] = fitcircle2(circleData)
  %[c0, r] = fitCircleRobust(circleData)
  %[c0, a, b] = fitellipseOnAxis2(circleData);
  %ERR_RATIO_XY =50/15;
  %%[c0, r] = fitCircleParametricRobust(circleData,ERR_RATIO_XY)
  %%plot the image for visual inspection
  %[srIm,yEdge_YZ,zEdge] =stormHist2d(y_rad{ii},z{ii},pixSize,'XLim',radLim,'YLim',radLim,'BlurSigma',sigmaBlur,'AmpScaleFactor',ampScaleFactor,'Gamma',gammaVal,'Normalize',satVal);
  %[srImXY,xEdge,yEdge_XY] =stormHist2d(x_th{ii},y_rad{ii},pixSize,'BlurSigma',sigmaBlur,'AmpScaleFactor',ampScaleFactor,'FlipY','Gamma',gammaVal,'Normalize',satVal);
  %subplot(2,1,1)
  %hold off
  %imagesc(yEdge_YZ,zEdge,srIm);
  %hold all;
  %%plotellipse(c0,r,r,0,'k-');
  %plotellipse(c0,a,b,0,'k-');

  %axis equal
  %subplot(2,1,2);
  %hold off
  %plot(meshRot{ii}(:,1),meshRot{ii}(:,2),'k-');
  %hold all;
  %plot(meshRot{ii}(:,3),meshRot{ii}(:,4),'k-');
  %imagesc(xEdge,yEdge_XY,srImXY);
  %cLine =  [mean(meshRot{ii}(:,[1,3]),2),mean(meshRot{ii}(:,[2,4]),2)];
  %plot(cLine(:,1),cLine(:,2),'k-');
  %axis equal;
  %pause
end
  

