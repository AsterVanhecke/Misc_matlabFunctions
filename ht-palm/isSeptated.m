function cell0 = isSeptated(cell0,sliceWidth,sliceStep,lenRange,plotOn)

mesh = cell0.mesh;
xcol = findSRField(cell0.localizations.parInfo,'semantic','position in sample space in X'); 
ycol = findSRField(cell0.localizations.parInfo,'semantic','position in sample space in Y'); 
x = cell0.localizations.data(:,xcol);
y = cell0.localizations.data(:,ycol);
minD= cell0.diameter.minD;
l0= minD(1);

[sliceN,sliceX]=getAlongCentreSlice(x,y,mesh,l0,sliceWidth,sliceStep,lenRange,plotOn);
cell0.septumPlot.sliceN = sliceN;
cell0.septumPlot.sliceX = sliceX;

%-----------------------------------------
function [sliceN,sliceX]=getAlongCentreSlice(x,y,mesh,l0,sliceWidth,sliceStep,lenRange,plotOn)
nBreak= 10;
LSTEP=50;

%recaluclate the fitting centreline
cLine = zeros(size(mesh,1),2);
cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);
[xC,yC]=spline2d(cLine(:,1),cLine(:,2),nBreak);
cLine=[xC(:),yC(:)];
len = sqrt((cLine(:,1)-cLine(1,1)).^2+(cLine(:,2)-cLine(1,2)).^2);% distances from the beginning of the mesh

%calculate the x,y position of the division site
cPt(:,1) = interp1(len,cLine(:,1),[l0-LSTEP/2, l0, l0+LSTEP/2]');
cPt(:,2) = interp1(len,cLine(:,2),[l0-LSTEP/2, l0, l0+LSTEP/2]');
c0 = cPt(2,:);%this is the midpoint

%get angle
dX = cPt(3,1)-cPt(1,1);
dY = cPt(3,2)-cPt(1,2);
theta = atan2(dY,dX);

%centre the xy data
x = x-c0(1);
y = y-c0(2);
%rotate parallel to xAxis 
thetaRot = -theta;
R = [cos(thetaRot),-sin(thetaRot);...
     sin(thetaRot),cos(thetaRot)];
[XYrot] = (R*[x';y'])';
xRot = XYrot(:,1);
yRot = XYrot(:,2);

%limits
xLim=[-lenRange/2, lenRange/2];
yLim = [-sliceWidth/2, sliceWidth/2];
isOk = xRot>=xLim(1) & xRot<=xLim(2) & yRot>=yLim(1) & yRot<=yLim(2);
xRotCrop = xRot(isOk);
yRotCrop = yRot(isOk);
slicePos = xRotCrop;
[sliceN,sliceX]=hist(slicePos,xLim(1):sliceStep:xLim(2));

if plotOn
  subplot(2,1,1)
  hold off
  plot(xRot,yRot,'rx');
  hold all;
  plot(xRotCrop,yRotCrop,'bx');
  axis equal
  subplot(2,1,2)
  plot(sliceX,sliceN);
  pause
end

%-------------------------------------------------------------
function [xOut,yOut]=spline2d(x,y,nBreak);
nPt = numel(x);
nPtOut = 10*nPt;% upsample the output, but keep discrete

t=linspace(0,1,nPt);
tOut=linspace(0,1,nPtOut);
ppX = splinefit(t,x,nBreak);
ppY = splinefit(t,y,nBreak);
xOut = ppval(ppX,tOut);
yOut = ppval(ppY,tOut);

