function plotHtpalmCrossSection(srInMesh, movieName,lRange,gammaVal,satVal,zLim,pixSize);

%pixSize=20;sigma=10;

sigma=10;
xMin = -4000;xMax = 4000;
yMin = -400;yMax = 400;
zMin = zLim(1);zMax = zLim(2);
ampScaleFactor = 100;
m=xMin:pixSize:xMax;
n=yMin:pixSize:yMax;
p = zMin:pixSize:zMax;

cellLength = getCellLength(srInMesh);

nFullCell = 0;
nCell = numel(srInMesh);
cellEdge = zeros(nCell,2);
[dum, order]=sort(cellLength(:,1));

%box = [xMin, yMin, (xMax-xMin), (yMax-yMin)]; 
%maxDensity = getCellMaxDensityCrossSection(srInMesh,pixSize,lRange,n,p,yMin,yMax);
%normVal = maxDensity*.5;

for kk = 1:nCell
  length = cellLength(order(kk),1);
  ii = order(kk);
  %plot the cell
  srData = srInMesh{ii}.localizations;
  %box    = srInMesh{ii}.box;
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
  nFullCell = nFullCell +1;

  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

  %exclude points outside of the lRange
  inCell = (x>lRange(1) & x<lRange(2));
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);

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

  [srIm] =stormHist2d(y,z,pixSize,'XLim',[yMin,yMax],'YLim',[zMin,zMax],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor);
  fname = [movieName];
  %[srIm] =stormHist2d(x,y,pixSize,'XLim',[xMin,xMax],'YLim',[yMin,yMax],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor);
  %fname = ['testSlicing.tif'];
  kk
  if kk ==1 
     writeMode = 'overwrite';
  else
     writeMode = 'append';
  end
  imwrite(uint16(srIm),fname,'WriteMode',writeMode);
  %pause
end
%%nFullCell

%----------------------------------------------------------
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

