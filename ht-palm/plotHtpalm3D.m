function  [VOL_blur, hF,isoVal,xEdge, yEdge, zEdge]= plotHtpalm3D(srInMesh, cellNo,kymographPixSize,lRange,findSpot_xLim,gammaVal,satVal,zLim,pixSize,sigmaXY,sigmaZ,pixSize3d,varargin);

xMin = -4000;
xMax = 4000;
yMin = -400;
yMax = 400;
zMin = zLim(1);
zMax = zLim(2);
zScale = 1;
doSortCells=false;
ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'SortCells')
    doSortCells=true;
    ii = ii + 1;
  elseif strcmp(varargin{ii},'ZScale')
    zScale = varargin{ii+1};
    ii = ii+2;
  elseif strcmp(varargin{ii},'XLim')
    xLim = varargin{ii+1};
    xMin = xLim(1);
    xMax = xLim(2);
    ii = ii + 2;
  elseif strcmp(varargin{ii},'YLim')
    yLim = varargin{ii+1};
    yMin = yLim(1);
    yMax = yLim(2);
    ii = ii + 2;
  else
    ii = ii + 1;
  end
end

m =xMin:pixSize:xMax;
n =yMin:pixSize:yMax;
p = zMin:pixSize:zMax;

kymographEdges = xMin:kymographPixSize:xMax;

cellLength = getCellLength(srInMesh);

nCell = numel(srInMesh);

if doSortCells
   [dum, order]=sort(cellLength(:,1));
else
   order = 1:nCell;
end

ii = order(cellNo);
%plot the cell
srData = srInMesh{ii}.localizations;
%box    = srInMesh{ii}.box;
mesh   = srInMesh{ii}.meshAligned;
%only show points inside ymin rand and x slice lRange
x=srInMesh{ii}.l;
y=srInMesh{ii}.d;
xcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in x dimension');
ycol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in y dimension');
zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in z dimension');
xExt=srInMesh{ii}.localizations.data(:,xcol);
yExt=srInMesh{ii}.localizations.data(:,ycol);
z=srInMesh{ii}.localizations.data(:,zcol);
z = z*zScale;%RESCALE Z

inCell = (y>yMin & y<yMax);
x=x(inCell);
y=y(inCell);
z=z(inCell);
xExt=xExt(inCell);
yExt=yExt(inCell);

%recentre and flip the data along CofM
[x,y,mesh] = alignCM(x,y,mesh);

ringRange = getBrightRegion(x,kymographEdges,kymographPixSize,lRange,findSpot_xLim);
 
%exclude points outside of the lRange
inCell = (x>ringRange(1) & x<ringRange(2));
x=x(inCell);
y=y(inCell);
z=z(inCell);
xExt=xExt(inCell);
yExt=yExt(inCell);

%flip yaxis (which has Z on it) so z is up.
hF = gcf;
hold off;
[VOL_blur,hP,isoVal,xEdge, yEdge, zEdge]= Render_3DVol_simple(xExt,yExt,z,'ZLim',[zMin,zMax],'SigmaXY',sigmaXY,'SigmaZ',sigmaZ,'PixSize',pixSize3d);

%----------------------------------------------
function setFigDefaults()


hF =gcf;
hAx = gca;
set(hAx,'TickDir','out')
set(hAx,'Box','off');
set(hF,'Color','w')
set(hAx,'FontSize',14);
%set(hAx,'YTick',[])



%-----------------------------------------------------------------------------------------------
function [density] =renderSingleCell(XPosition,YPosition,n,m,pixSize,normVal)

minX = min(m);
maxX= max(m);
minY = min(n);
maxY= max(n);
%remove out of bounds data
isInBounds = XPosition > minX & XPosition < maxX ...
              & YPosition > minY & YPosition < maxY ;

XPosition = XPosition(isInBounds);
YPosition = YPosition(isInBounds);

subset = logical(ones(size(XPosition)));

nx = numel(m);
ny = numel(n);
notrials = 1; % SET NOTRIALS TO 20 FOR JITTERING
NN = numel(XPosition);
%amount = (sx+sy)/2;
amount = XPosition*0; % TO AVOID JITTERING
density = zeros(ny,nx);
pxx=pixSize; pxy=pixSize; 
pxVol=pxx * pxy ;

density = hist3([YPosition,XPosition],'Edges',{n,m});

density(density<0)=0;
%density = density/normVal;

%---------------------------------------------------------------------------
function maxDensity = getCellMaxDensityCrossSection(srInMesh,pixSize,lRange,n,p,yMin,yMax);

maxDensity = 0;
cellLength = getCellLength(srInMesh);
nCell = numel(srInMesh);
cellEdge = zeros(nCell,2);
[dum, order]=sort(cellLength(:,1));
for kk = 1:nCell
  length = cellLength(order(kk),1);
  ii = order(kk);
  %plot the cell
  srData = srInMesh{ii}.localizations;
  mesh   = srInMesh{ii}.meshAligned;
  %only show points inside ymin rand and x slice lRange
  x=srInMesh{ii}.l;
  y=srInMesh{ii}.d;
  zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in z dimension');
  z=srInMesh{ii}.localizations.data(:,zcol);

  inCell = (y>yMin & y<yMax);
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);

  %plot the cell
  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

  %exclude points outside of the lRange
  inCell = (x>lRange(1) & x<lRange(2));

  %make a box
  normVal =1;
  [srIm] =renderSingleCell(y,z,n,p,pixSize,normVal);
  imMaxDensity = max(srIm(:));
  if imMaxDensity >maxDensity
     maxDensity = imMaxDensity
  end
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

%hold off;
%plot(kymographEdges,nK,'k');
%hold all;
%plot([xCM,xCM],[0,nKMax],'r');
%pause;
