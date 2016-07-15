function [ringRange, nMolecules, locsRingDZ] =  plotHtpalmCrossSection2(srInMesh, movieName,kymographPixSize,lRange,findSpot_xLim,gammaVal,satVal,zLim,pixSize,sigma,varargin);

xMin = -4000;xMax = 4000;
yMin = -400;yMax = 400;
zMin = zLim(1);zMax = zLim(2);
ampScaleFactor = 100;
doSortCells=false;
doNorm = true;
zScale = 1;
ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'SortCells')
    doSortCells=true;
    ii = ii + 1;
  elseif strcmp(varargin{ii},'ZScale')
    zScale = varargin{ii+1};
    ii = ii+2;
  elseif strcmp(varargin{ii},'Normalize')
    doNorm = varargin{ii+1};
    ii = ii+2;
  else
    ii = ii + 1;
  end
end

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

for kk = 1:nCell
  ii = order(kk);
  %plot the cell
  srData = srInMesh{ii}.localizations;
  %box    = srInMesh{ii}.box;
  mesh   = srInMesh{ii}.meshAligned;
  %only show points inside ymin rand and x slice lRange
  x=srInMesh{ii}.l;
  y=srInMesh{ii}.d;
  zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in Z');
  z=srInMesh{ii}.localizations.data(:,zcol);
  %rescale z
  z = z*zScale;

  inCell = (y>yMin & y<yMax);
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);

  %plot the cell

  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

  kk
  ringRange(kk,:) = getBrightRegion(x,kymographEdges,kymographPixSize,lRange,findSpot_xLim);
   
  %exclude points outside of the lRange
  inCell = (x>ringRange(kk,1) & x<ringRange(kk,2));
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);

  locsRingDZ{kk} = [y(:),z(:)];
  nMolecules(kk) = numel(y);

  %make a box
  %extraSpace = 100;%nm
  %xMax = max([mesh(:,1);mesh(:,3)])+extraSpace;
  %yMin = min([mesh(:,2);mesh(:,4)])-extraSpace;
  %yMax = max([mesh(:,2);mesh(:,4)])+extraSpace;
  %[srIm] =renderSingleCell(y,z,n,p,pixSize,normVal);
  %hold off
  %imagesc(n,p,srIm,[0 1]);
  %colormap(hot);
  %hold all

  %axis equal
  %ylim([zMin, zMax]);
  %xlim([yMin, yMax]);
  %setFigDefaults();
  %plot(srData.data(:,1),srData.data(:,2),'rx');
  %box
  %h1=gcf;

  %flip yaxis (which has Z on it) so z is up.
  if doNorm 
    [srIm] =stormHist2d(y,z,pixSize,'XLim',[yMin,yMax],'YLim',[zMin,zMax],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor,'FlipY','Gamma',gammaVal,'Normalize',satVal);
  else
    [srIm] =stormHist2d(y,z,pixSize,'XLim',[yMin,yMax],'YLim',[zMin,zMax],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor,'FlipY','Gamma',gammaVal);
  end
  
  fname = [movieName];
  if kk ==1 
     writeMode = 'overwrite';
  else
     writeMode = 'append';
  end

  if doNorm 
    %normalized images are between 0->1 so need to stretch to uint16
    imwrite(im2uint16(srIm),fname,'WriteMode',writeMode);
  else %raw images are already integer counts so just round and cast
    imwrite(cast(round(srIm),'uint16'),fname,'WriteMode',writeMode);
  end
end
save([fname(1:end-3),'_ringRange.mat'],'ringRange');
%%nFullCell

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
  zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in Z');
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
