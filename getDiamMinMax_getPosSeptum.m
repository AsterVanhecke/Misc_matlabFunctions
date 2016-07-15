function [maxD minD] = getDiamMinMax(cell0,plotOn)

if ~exist('plotOn','var')
  plotOn=false;
end
nBreak = 10;
len = cell0.diameter.length;
diam_fwhm = cell0.diameter.diam_fwhm;
diam_mesh = cell0.diameter.diam_mesh;

%to find the max diameter - needs to be robust so use spline smoothing
%nPt = numel(len);
ppDiam = splinefit(len,diam_fwhm,nBreak);
diamSmth=ppval(ppDiam,len);

%find the first local maximum on the left
[posLeftStart,valLeftStart,idxLeftStart]   =  getLeft(len,diamSmth);
[posRightEnd,valRightEnd,idxRightEnd] = getRight(len,diamSmth);
if isempty(posLeftStart)|| isempty(posRightEnd)
  maxD=[];
  minD=[];
else
  % then on the right.
  %this marks the limits of the search area
  %find the absolute minimum
  [minVal idxMin ] = min(diam_fwhm(idxLeftStart:idxRightEnd));
  idxMin = idxMin+idxLeftStart-1;
  minD = [len(idxMin),minVal, idxMin];
  %find the absolute maximum to the the left and right of the minima
  %left
  [maxVal idxMax] = max(diamSmth(idxLeftStart:idxMin));
  idxMax= idxMax + idxLeftStart-1;
  maxD(1,:) = [len(idxMax),maxVal];
  %Right
  [maxVal idxMax] = max(diamSmth(idxMin:idxRightEnd));
  idxMax= idxMax +idxMin -1;
  maxD(2,:) = [len(idxMax),maxVal];


  if plotOn
    hold off;
    plot(len,diam_fwhm);
    hold all;
    plot(len,diamSmth);
    plot(posLeftStart,valLeftStart,'ro');
    plot(posRightEnd,valRightEnd,'ro');
    plot(minD(:,1),minD(:,2),'ko');
    plot(maxD(:,1),maxD(:,2),'ko');
  end
end

%---------------------------------------------------------
function [posLeftMax,valLeftMax,idxLeftMax] = getLeft(x,y);
nX = numel(x);
ii =2;
foundLeftPk = false;
idxLeftMax=[];
valLeftMax=[];
posLeftMax=[];
while ~foundLeftPk && ii<=nX-1
  if (y(ii)-y(ii-1)>0)&&(y(ii+1)-y(ii)<0)%absolute maximum
    idxLeftMax = ii;
    valLeftMax = y(ii);
    posLeftMax=x(ii);
    foundLeftPk=true;
 % elseif (y(ii)-y(ii-1)<0)%if we're descending we're past the peak
 %   idxLeftMax = 1;
 %   valLeftMax = y(1);
 %   posLeftMax=x(1);
 %   foundLeftPk=true;
  else
    ii=ii+1;
  end
end

%---------------------------------------------------------
function [posRightMax,valRightMax,idxRightMax] = getRight(x,y);
nX = numel(x);
ii =nX-1;
foundRightPk = false;
idxRightMax=[];
valRightMax=[];
posRightMax=[];
while ~foundRightPk && ii>=2
  if (y(ii)-y(ii-1)>0)&&(y(ii+1)-y(ii)<0)
    idxRightMax = ii;
    valRightMax = y(ii);
    posRightMax=x(ii);
    foundRightPk=true;
%  elseif (y(ii)-y(ii-1)>0)%if we're descending we're past the peak
%    idxRightMax =nX;
%    valRightMax = y(nX);
%    posRightMax=x(nX);
%    foundRightPk=true;
  else
    ii=ii-1;
  end
end


