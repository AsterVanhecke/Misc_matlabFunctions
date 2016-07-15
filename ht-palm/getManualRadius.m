function [radius circleInfo isGoodCell] = getManualRadius(srInMesh,kymographPixSize,lRange,findSpot_xLim,gammaVal,satVal,zLim,pixSize,sigmaBlur)

close all;

xMin = -4000;xMax = 4000;
yMin = -400;yMax = 400;
zMin = zLim(1);zMax = zLim(2);
ampScaleFactor = 100;
doSortCells=false;
zScale = 1;
ii = 1;

m=xMin:pixSize:xMax;
n=yMin:pixSize:yMax;
p = zMin:pixSize:zMax;
kymographEdges = xMin:kymographPixSize:xMax;

cellLength = getCellLength(srInMesh);

nCell = numel(srInMesh);

if doSortCells
   [dum, order]=sort(cellLength(:,1));
else
   order = 1:nCell;
end

isGoodCell = zeros(nCell,1);
h1=figure;
h2 = figure;
for kk = 1:nCell
  ii = order(kk);
  ii
  %plot the cell
  srData = srInMesh{ii}.localizations;
  %box    = srInMesh{ii}.box;
  mesh   = srInMesh{ii}.meshAligned;
  %only show points inside ymin rand and x slice lRange
  x=srInMesh{ii}.l;
  y=srInMesh{ii}.d;
  l0=x;d0=y;
  zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in Z');
  z=srInMesh{ii}.localizations.data(:,zcol);
  %rescale z
  z = z*zScale;

  inCell = find(y>yMin & y<yMax);
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);

  %plot the cell

  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

  kk
  ringRange(kk,:) = getBrightRegion(x,kymographEdges,kymographPixSize,lRange,findSpot_xLim);
   
  %exclude points outside of the lRange
  inRing= find(x>ringRange(kk,1) & x<ringRange(kk,2));
  x0=x;y0=y;
  x=x(inRing);
  y=y(inRing);
  z=z(inRing);
  inRingAllLoc = inCell(inRing);%this list of coords works for raw xyz localizations (ie not on the shortened list)

  locsRingDZ{kk} = [y(:),z(:)];
  nMolecules(kk) = numel(y);

  figure(h1);
  extraSpace = 100;%nm
  normVal =1;
  xMax = max([mesh(:,1);mesh(:,3)])+extraSpace;
  yMin = min([mesh(:,2);mesh(:,4)])-extraSpace;
  yMax = max([mesh(:,2);mesh(:,4)])+extraSpace;
  plotSingleCell(srData, srInMesh{ii}.mesh,inRingAllLoc,x0,y0,x,y);

  %flip yaxis (which has Z on it) so z is up.
  [srIm,d,z] =stormHist2d(y,z,pixSize,'XLim',[yMin,yMax],'YLim',[zMin,zMax],'BlurSigma',sigmaBlur,'AmpScaleFactor',ampScaleFactor,'FlipY','Gamma',gammaVal,'Normalize',satVal);
  
  figure(h2);
  imagesc(d,z,srIm);
  colormap(hot);
  set(gca,'YDir','normal');
  axis equal

  inputOk = false;
  while ~inputOk
    str = input('Good cell? (''g''/''b''):','s');
    if strcmp(str,'g')
      isGoodCell(ii)  = 1;
      [radius{ii}, circleInfo{ii}] =getCircle(h2);
      inputOk = true;
    elseif strcmp(str,'b')
      isGoodCell(ii)  = 0;
      radius{ii}=[];
      circleInfo{ii}=[];
      inputOk = true;
    end
  end
  save('tmpRingResult150215_1.mat','isGoodCell','radius','circleInfo');
end

%---------------------------------------------------------
function [xRange, xCM] = getBrightRegion(x,kymographEdges,kymographPixSize,lRange,findSpot_xLim);
x = x(x>findSpot_xLim(1) & x<findSpot_xLim(2));%trip to within allowed region

%first, plot a 1d intensity plot (along x)
nK = hist(x,kymographEdges);
[nKMax,idx] = max(nK);
xNMax = kymographEdges(idx);

%refine position of bright spot by calculating CM in range lRange;
inXRegion = find( kymographEdges>=xNMax-lRange & kymographEdges<=xNMax+lRange);
xInReg = kymographEdges(inXRegion);
nInReg = nK(inXRegion);

xCM = sum(xInReg.*nInReg)/sum(nInReg);
xRange = [(xCM -lRange),(xCM +lRange)];

%-----------------------------------------------------------------------------------------------
function plotSingleCell(srData,mesh,inRing,l0,d0,l,d)

xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');

x = srData.data(:,xcol);
y = srData.data(:,ycol);
xRing = x(inRing);
yRing = y(inRing);

subplot(2,1,1)
hold off;
plot(mesh(:,1),mesh(:,2),'k-')%left side of cell
hold all;
plot(mesh(:,3),mesh(:,4),'k-');%right side of cell
plot(x,y,'r.');
plot(xRing,yRing,'g.');
axis equal
subplot(2,1,2)
hold off;
plot(l0,d0,'r.')
hold all;
plot(l,d,'g.')
axis equal

%-----------------------------------------------------------------------------------------------
function [radius, circleCentre] =getCircle(h2)
hC = imellipse(get(h2,'Children'));
str = input('Press Return when finished selecting ring.','s');

pos = getPosition(hC);
radius = pos(3:4)
circleCentre = [pos(1)+radius(1)/2,pos(1)+radius(1)/2 ]
