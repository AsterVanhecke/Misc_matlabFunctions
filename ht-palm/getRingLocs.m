function [x_th,y_rad,z_ring,t, meshRot] = getRingLocs(srInMesh,kymographPixSize,lRange,findSpot_xLim,gammaVal,satVal,zLim,pixSize,sigmaBlur)

lMin = -4000;lMax = 4000;
dMin = -400;dMax = 400;
zMin = zLim(1);zMax = zLim(2);
ampScaleFactor = 100;
doSortCells=false;
zScale = 1;
ii = 1;

m=lMin:pixSize:lMax;
n=dMin:pixSize:dMax;
p = zMin:pixSize:zMax;
kymographEdges = lMin:kymographPixSize:lMax;

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
  l=srInMesh{ii}.l;
  d=srInMesh{ii}.d;
  xcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in X');
  ycol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in Y');
  zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in Z');
  x=srInMesh{ii}.localizations.data(:,xcol);
  y=srInMesh{ii}.localizations.data(:,ycol);
  z=srInMesh{ii}.localizations.data(:,zcol);
  %rescale z
  z = z*zScale;

  %exclude points outside of the width range
  inCell = find(d>dMin & d<dMax);
  l=l(inCell);
  d=d(inCell);
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);
  l0=l;d0=d;
  %exclude points outside of the Z-range
  inCell = find(z>zLim(1)& z<zLim(2));
  l=l(inCell);
  d=d(inCell);
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);

  %plot the cell

  %recentre and flip the data along CofM
  [l,d,mesh] = alignCM(l,d,mesh);

  kk
  ringRange(kk,:) = getBrightRegion(l,kymographEdges,kymographPixSize,lRange,findSpot_xLim);
   
  %exclude points outside of the lRange
  inRing= find(l>ringRange(kk,1) & l<ringRange(kk,2));
  l=l(inRing);
  d=d(inRing);
  x=x(inRing);
  y=y(inRing);
  z=z(inRing);
  inRingAllLoc = inCell(inRing);%this list of coords works for raw xyz localizations (ie not on the shortened list)

  locsRingDZ{kk} = [d(:),z(:)];
  nMolecules(kk) = numel(d);

  steplength = srInMesh{kk}.steplength*pixSize;
  [x_th{kk} y_rad{kk} meshRot{kk}] = getRotatedRing(ringRange(kk,:),x,y,mesh,pixSize,steplength,srInMesh{kk});

  t(kk) =srInMesh{kk}.tExpt_min; 
  z_ring{kk} = z;

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
