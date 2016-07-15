function [hFig srIm locs nMol] = plotHtpalm2D(srInMesh,cellNo,pixSize,sigma,varargin);

extraSpace = 100;%nm
normVal = 1000;
ampScaleFactor = 100;
doSortCells=false;
autoXLim = true;
autoYLim = true;
tZero = 0;%offset to account for time after start of synchrony
includeTime=true;
stormHistArg = {};
ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'SortCells')
    doSortCells=true;
   includeTime=false;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'NormVal')
   normVal = varargin{ii+1};
    ii = ii + 2;
 elseif strcmp(varargin{ii},'Normalize')
   stormHistArg = {stormHistArg{:},varargin{ii:ii+1}};
   ii = ii + 2;
 elseif strcmp(varargin{ii},'Gamma')
   stormHistArg = {stormHistArg{:},varargin{ii:ii+1}};
   ii = ii + 2;
 elseif strcmp(varargin{ii},'XLim')
    xLim = varargin{ii+1};
    minX = xLim(1);
    maxX = xLim(2);
    autoXLim=false;
    ii = ii + 2;
 elseif strcmp(varargin{ii},'YLim')
    yLim = varargin{ii+1};
    minY = yLim(1);
    maxY = yLim(2);
    autoYLim=false;
    ii = ii + 2;
  elseif strcmp(varargin{ii},'TStart')
    tZero= varargin{ii+1};
    ii = ii + 2;
 else 
   ii=ii+1;
 end
end

nCell = numel(srInMesh);
cellLength = getCellLength(srInMesh);
if doSortCells
   [dum, order]=sort(cellLength(:,1));
else
   order = 1:nCell;
end

ii = order(cellNo)
%plot the cell
srData = srInMesh{ii}.localizations;
%box    = srInMesh{ii}.box;
mesh   = srInMesh{ii}.mesh;
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
x = srData.data(:,xcol);
y = srData.data(:,ycol);

%recentre and flip the data along CofM
[xC,yC,meshC] = alignHorizontal(x,y,mesh);
[xC,yC,minX_tmp,maxX_tmp,minY_tmp,maxY_tmp] = clipData(xC,yC,meshC,extraSpace);

if autoXLim
  minX = minX_tmp;
  maxX = maxX_tmp;
end

if autoYLim
  minY = minY_tmp;
  maxY = maxY_tmp;
end

%save the non-normalised images
[srIm,xEdge,yEdge] =stormHist2d(xC,yC,pixSize,'XLim',[minX,maxX],'YLim',[minY,maxY],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor,stormHistArg{:});

hold off
imagesc(xEdge,yEdge,srIm,[0 1]);
set(gca,'YDir','normal');
colormap(hot);
hold all
plot(meshC(:,1),meshC(:,2),'y-')%left side of cell
plot(meshC(:,3),meshC(:,4),'y-');%right side of cell
axis equal
ylim([minY,maxY]);
xlim([minX,maxX]);
if includeTime
   t = srInMesh{ii}.tExpt_min
   title(['# ',num2str(cellNo),', ',sprintf('%.1f',t+tZero),'min']);
else
   title(['# ',num2str(cellNo)]);
end
setFigDefaults();
hFig=gcf;
%fname = [figFolder,filesep,figName,'_', sprintf('%04d',kk),'.tif'];
%saveas(h1,fname);
locs = [xC(:),yC(:)];
nMol = numel(xC);
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

%-------------------------------------------------------
function [x0,y0,minX,maxX,minY,maxY] = clipData(x,y,mesh,extraSpace);

minX = min([mesh(:,1);mesh(:,3)]) - extraSpace;
maxX = max([mesh(:,1);mesh(:,3)]) + extraSpace;
minY = min([mesh(:,2);mesh(:,4)])- extraSpace;
maxY = max([mesh(:,2);mesh(:,4)])+ extraSpace;

locIsOk = (x>=minX & x<=maxX & y>=minY & y<=maxY);
x0 = x(locIsOk);
y0 = y(locIsOk);

%----------------------------------
