function [width, xLeft, xRight,yLeft,yRight]= getFWHM(x,y,xLim,searchStep)
%find first maximum on left

[posLeft,posLeftMax,valLeftMax] = getLeft(x,y,xLim(1),searchStep);
[posRight,posRightMax,valRightMax] = getRight(x,y,xLim(2),searchStep);

%%DEBUG
%hold off;
% plot(x,y);
% hold all;
% plot(posLeftMax,valLeftMax,'ro');
% plot(posLeft,valLeftMax/2,'bo')
% plot(posRightMax,valRightMax,'ro');
% plot(posRight,valRightMax/2,'bo');
%%DEBUG

width = posRight-posLeft;
xLeft = posLeft;
xRight = posRight;
yLeft =valLeftMax/2;
yRight =valRightMax/2;

%---------------------------------------------------------
function [posLeft,posLeftMax,valLeftMax] = getLeft(x,y,xLim,searchStep)
nX = numel(x);
ii =2;
foundLeftPk = false;
idxLeftMax=[];
valLeftMax=[];
posLeftMax=[];
while ~foundLeftPk && ii<=nX-1
  if (y(ii)-y(ii-1)>0)&&(y(ii+1)-y(ii)<=0)
    idxLeftMax = ii;
    valLeftMax = y(ii);
    posLeftMax=x(ii);
    foundLeftPk=true;
  elseif (y(ii)-y(ii-1)<0)%if we're descending we're past the peak
    idxLeftMax = 1;
    valLeftMax = y(1);
    posLeftMax=x(1);
    foundLeftPk=true;
  else
    ii=ii+1;
  end
end

if isempty(idxLeftMax)
  posLeft = NaN;
  posLeftMax=NaN;
  valLeftMax=NaN;
else
  %find the half max value
  posLeft = x(idxLeftMax);
  valCur = Inf;
  while posLeft>=x(1) && valCur > valLeftMax/2
    posLeft =posLeft-searchStep;
    valCur = interp1(x,y,posLeft);
  end
end

%---------------------------------------------------------
function [posRight,posRightMax,valRightMax] = getRight(x,y,xLim,searchStep)
nX = numel(x);
ii =nX-1;
foundRightPk = false;
idxRightMax=[];
valRightMax=[];
posRightMax=[];
while ~foundRightPk && ii>=2
  if (y(ii)-y(ii-1)>=0)&&(y(ii+1)-y(ii)<0)
    idxRightMax = ii;
    valRightMax = y(ii);
    posRightMax=x(ii);
    foundRightPk=true;
  elseif (y(ii)-y(ii-1)>0)%if we're descending we're past the peak
    idxRightMax =nX;
    valRightMax = y(nX);
    posRightMax=x(nX);
    foundRightPk=true;
  else
    ii=ii-1;
  end
end

if isempty(idxRightMax)
  posRight = NaN;
  posRightMax=NaN;
  valRightMax=NaN;
else
  %find the half max value
  posRight = x(idxRightMax);
  valCur = Inf;
  while posRight<=x(end) && valCur > valRightMax/2
    posRight =posRight+searchStep;
    valCur = interp1(x,y,posRight);
  end
end

