function plotHtpalm2DInt(srInMesh,movieName,gammaVal,satVal,figFolder,figName,pixSize,sigma,varargin);

saveRawIm = pause;
saveMatFig = true;
normVal = 1000;
%pixSize=20;sigma=10;
xMin = -4000;xMax = 4000;
yMin = -400;yMax = 400;
ampScaleFactor = 100;
doSortCells=false;
ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'SortCells')
    doSortCells=true;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'NormVal')
   normVal = varargin{ii+1};
    ii = ii + 2;
  else
    ii = ii + 1;
  end
end

n=xMin:pixSize:xMax;
m=yMin:pixSize:yMax;


nFullCell = 0;
nCell = numel(srInMesh);

box = [xMin, yMin, (xMax-xMin), (yMax-yMin)]; 

cellLength = getCellLength(srInMesh);
if doSortCells
   [dum, order]=sort(cellLength(:,1));
else
   order = 1:nCell;
end

for kk = 1:nCell
  ii = order(kk)
  t = srInMesh{ii}.tExpt_min
  %plot the cell
  srData = srInMesh{ii}.localizations;
  %box    = srInMesh{ii}.box;
  mesh   = srInMesh{ii}.meshAligned;
  %only show points inside ymin rand and x slice lRange
  x=srInMesh{ii}.l;
  y=srInMesh{ii}.d;
  inCell = (y>yMin & y<yMax);
  x=x(inCell);
  y=y(inCell);

  %plot the cell
  nFullCell = nFullCell +1;

  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

  %make a box
  extraSpace = 100;%nm
  %xMin = min([mesh(:,1);mesh(:,3)])-extraSpace;
  %xMax = max([mesh(:,1);mesh(:,3)])+extraSpace;
  %yMin = min([mesh(:,2);mesh(:,4)])-extraSpace;
  %yMax = max([mesh(:,2);mesh(:,4)])+extraSpace;
  %[srIm] =renderSingleCell(x,y,m,n,pixSize,normVal);
  %hold off
  %imagesc(n,m,srIm,[0 1]);
  %colormap(hot);
  %hold all
  %plot(mesh(:,1),mesh(:,2),'y-')%left side of cell
  %plot(mesh(:,3),mesh(:,4),'y-');%right side of cell
  %plot(mesh(:,3),mesh(:,4),'y-');%right side of cell

  %axis equal
  %ylim([yMin, yMax]);
  %xlim([xMin, xMax]);
  %setFigDefaults();
  %%plot(srData.data(:,1),srData.data(:,2),'rx');
  %%box
  %h1=gcf;
  %saveas(h1,fname);
  %%pause

  %save the non-normalised images
  [srIm,xEdge,yEdge] =stormHist2d(x,y,pixSize,'XLim',[xMin,xMax],'YLim',[yMin,yMax],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor);
  kk
  if saveRawIm == true
     if kk ==1 
        writeMode = 'overwrite';
     else
        writeMode = 'append';
     end
     imwrite(uint16(srIm),movieName,'WriteMode',writeMode);
  end
  if saveMatFig ==true
     if kk ==1 & ~exist(figFolder,'dir')
        mkdir(figFolder);
     end

     srIm = srIm./normVal;
     srIm(srIm>1) = 1;
     hold off
     imagesc(xEdge,yEdge,srIm,[0 1]);
     colormap(hot);
     hold all
     plot(mesh(:,1),mesh(:,2),'y-')%left side of cell
     plot(mesh(:,3),mesh(:,4),'y-');%right side of cell
     plot(mesh(:,3),mesh(:,4),'y-');%right side of cell
     ylim([yMin, yMax]);
     xlim([xMin, xMax]);
     axis equal
     setFigDefaults();
     h1=gcf;
     fname = [figFolder,filesep,figName,'_', sprintf('%04d',kk),'.tif'];
     saveas(h1,fname);
  end
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
function [density] =renderSingleCell(XPosition,YPosition,m,n,pixSize,normVal)

minX = min(n);
maxX= max(n);
minY = min(m);
maxY= max(m);
%remove out of bounds data
isInBounds = XPosition > minX & XPosition < maxX ...
              & YPosition > minY & YPosition < maxY ;

XPosition = XPosition(isInBounds);
YPosition = YPosition(isInBounds);

subset = logical(ones(size(XPosition)));

nx = numel(n);
ny = numel(m);
notrials = 1; % SET NOTRIALS TO 20 FOR JITTERING
NN = numel(XPosition);
%amount = (sx+sy)/2;
amount = XPosition*0; % TO AVOID JITTERING
density = zeros(ny,nx);
pxx=pixSize; pxy=pixSize; 
pxVol=pxx * pxy ;

density = hist3([YPosition,XPosition],'Edges',{m,n});

density(density<0)=0;
density = density/normVal;

%---------------------------------------------------------------------------
function maxDensity = getCellMaxDensity(srInMesh,pixSize,n,m);

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
  %only show points inside cell
  IN = srInMesh{ii}.position>0;
  x=srInMesh{ii}.l(IN);
  y=srInMesh{ii}.d(IN);
  %plot the cell

  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

  %make a box
  normVal =1;
  [srIm] = renderSingleCell(x,y,m,n,pixSize,1);
  imMaxDensity = max(srIm(:));
  if imMaxDensity >maxDensity
     maxDensity = imMaxDensity
  end
end

