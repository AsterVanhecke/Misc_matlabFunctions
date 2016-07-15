function renderKymographLengthAvg(srInMesh,figName,lenStepSz,varargin)

doSortCells = false;
doNormLength = false;
pixSize=100;
sigma=10;
gammaVal = .3;
satVal=0;
is3D = false;
ii = 1;
sumAxis = 1;
scaleFac = 2000;
normIntensity = false;
xRange = [-2000 2000];
while ii <= numel(varargin)
  if strcmp(varargin{ii},'Pixel Size')
    pixSize = varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Sigma')
    sigma= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Gamma')
    gammaVal= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Saturate')
    satVal= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'NormaliseIntensity')
    normIntensity= true;
    ii = ii + 1;
  elseif strcmp(varargin{ii},'PlotX')
    sumAxis = 1;
    ii = ii + 1;
  elseif strcmp(varargin{ii},'PlotY')
    sumAxis = 2;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'SortCells')
    doSortCells=true;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'NormaliseLength')
    doNormLength = true
    ii = ii + 1;
  elseif strcmp(varargin{ii},'XRange')
    xRange= varargin{ii+1};
    ii = ii + 2;
  else
     ii=ii+1;
  end
end


%make a box
extraSpace = 100;%nm
%xMin = min([mesh(:,1);mesh(:,3)])-extraSpace;
%xMax = max([mesh(:,1);mesh(:,3)])+extraSpace;
%yMin = min([mesh(:,2);mesh(:,4)])-extraSpace;
%yMax = max([mesh(:,2);mesh(:,4)])+extraSpace;
if ~doNormLength
  xMin = xRange(1);
  xMax = xRange(2);
else
   xMin = -0.75*scaleFac;
   xMax = 0.75*scaleFac;
end

yMin = -400;yMax = 400;
box = [xMin, yMin, (xMax-xMin), (yMax-yMin)];

nCell = numel(srInMesh);
srIm = zeros(nCell,numel(xMin:pixSize:xMax));
cellEdge = zeros(nCell,2);

for kk = 1:nCell
  %plot the cell
  srData = srInMesh{ii}.localizations;
  %box    = srInMesh{kk}.box;
  mesh   = srInMesh{kk}.meshAligned;
  %only show points inside ymin
  x=srInMesh{kk}.l;
  y=srInMesh{kk}.d;
  inCell = (y>-400 & y<400);
  x=x(inCell);
  y=y(inCell);

  %recentre the data around the middle of the cell
  [x,y,mesh] = alignCM(x,y,mesh);
  %normalise length
  norm = @(x) (x - min(x))/(max(x)-min(x));
  if doNormLength
     x = norm(x)*scaleFac - scaleFac/2;
  end

  [imTmp,m,n] = hist2DKym(x,y,pixSize,sigma,sumAxis, box);
  srIm(kk,:) = imTmp;
   if doNormLength
     cellEdge(kk,:) = [-scaleFac/2 +scaleFac/2];
   else
     cellEdge(kk,:) = [min([mesh(:,1);mesh(:,3)]) max([mesh(:,1);mesh(:,3)]) ];
   end
end

length = getLengthExpt(srInMesh);

lenPcVal = 0:lenStepSz:1;
nStep = numel(lenPcVal);
for ii = 1:nStep
   lenVal(ii) = percentileLim(length,lenPcVal(ii));
end

imWidth = size(srIm,2);
srImLenAvg = zeros(nStep,imWidth);

for ii = 1:nStep-1
   lenOk = (length>= lenVal(ii) & length < lenVal(ii+1));
   srImLenAvg(ii,:) = mean(srIm(lenOk,:),1);
   cellEdgeLenAvg(ii,:) = mean(cellEdge(lenOk,:),1);
end

lenPlot = (lenVal(2:end) - (lenVal(2) - lenVal(1))/2 ) /1e3;

%adjust kymograph image
srImLenAvg= srImLenAvg/max(srImLenAvg(:));
srImLenAvg = saturateImage(srImLenAvg,satVal);
srImLenAvg = adjustGamma(srImLenAvg,gammaVal);

h1=figure;
xpos_um = n/1e3;
cellEdge_nm = cellEdgeLenAvg/1e3;
cellEdge_nm(:,1) = smooth(cellEdge_nm(:,1),10,'rloess');
cellEdge_nm(:,2) = smooth(cellEdge_nm(:,2),10,'rloess');
hIm = imagesc(xpos_um,lenPlot,srImLenAvg);

hold all;
%plot(cellEdge_nm(:,1),lenPlot,'k-','LineWidth',2);
%plot(cellEdge_nm(:,2),lenPlot,'k-','LineWidth',2);
ax1 = gca;

load('colormap_gg.mat');
colormap(colormap_gg);

%caxis([0, max(srImTimeAvg(:))]);

setFigDefaults();
set(ax1,'TickDir', 'out');
set(ax1,'YDir','Reverse')
%set(gca,'XTick',[-2.0 0.0 2.0])
%set(gca,'YTick',[0 100 200 300])
hXLabel = xlabel('Distance ({\mu}m)');
hYLabel = ylabel('Cell length ({\mu}um)');
set(gca, 'LineWidth', 2)
set(gcf,'PaperPositionMode','auto');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');

saveas(h1,figName);
%print('-depsc2',[figName(1:end-3),'eps']);

%-----------------------------------------------------------------------------------------------
function [density,m,n] =hist2DKym(XPosition,YPosition,pixSize,sigma,sumAxis,box)

if ~exist('box','var')
  minX = min(XPosition);
  maxX = max(XPosition);
  minY = min(YPosition);
  maxY = max(YPosition);
else
  minX = box(1);
  maxX = box(1)+box(3);
  minY = box(2);
  maxY = box(2)+box(4);
end

%remove out of bounds data
isInBounds = XPosition > minX & XPosition < maxX ...
              & YPosition > minY & YPosition < maxY ;

XPosition = XPosition(isInBounds);
YPosition = YPosition(isInBounds);

sx = ones(size(XPosition))*sigma;
sy=sx;
subset = logical(ones(size(XPosition)));

%n=linspace(minX,maxX,res);
%m=linspace(minY,maxY,res);
n=minX:pixSize:maxX;
m=minY:pixSize:maxY;
%n=maxX:-pixSize:minX;%somehow the axes end up getting flipped otherwise
%m=maxY:-pixSize:minY;%somehow the axes end up getting flipped otherwise
nx = numel(n);
ny = numel(m);
notrials = 1; % SET NOTRIALS TO 20 FOR JITTERING
NN = numel(XPosition);
amount = (sx+sy)/2;
amount = amount*0; % TO AVOID JITTERING
density = zeros(ny,nx);
pxx=pixSize; pxy=pixSize; 
pxVol=pxx * pxy ;
for i=1:notrials
%    RR = [round((res-1)*[...
%        (XPosition(subset)+randn(NN,1).*amount(subset)-minX)/(maxX-minX) ...
%        (YPosition(subset)+randn(NN,1).*amount(subset)-minY)/(maxY-minY) ])+1 ...
%          round(nz*(ZPosition(subset)+randn(NN,1).*amount(subset)-minZ)/(maxZ-minZ))];
%    RRok = all(RR(:,1:2)<=res,2) & all(RR>=1,2) & RR(:,3)<=nz;
    RR = [...
          round((ny-1)*(YPosition(subset)+randn(NN,1).*amount(subset)-minY)/(maxY-minY))+1 ...
          round((nx-1)*(XPosition(subset)+randn(NN,1).*amount(subset)-minX)/(maxX-minX))+1];
    d = accumarray(RR,1,[ny,nx])/pxVol;
    if i==1
        density = d;
    else
        density = density + d;
    end

    density = density/notrials;
    
end

density = sum(density,sumAxis);
%TODO use a faster gauss filter here
sPix = sigma/pxx;
gWindow = ceil(5*sPix);
gKern = fspecial('gaussian',gWindow, sPix);
dMax = max(density(:));
density = imfilter(density,gKern,'replicate');

density(density<0)=0;

%-----------------------------------------------------------------------------------------------
function imG= adjustGamma(im,gammaVal)

imMax = max(im(:));
%normalise image
imG = ((im/imMax).^gammaVal)*imMax;

%-----------------------------------
%function b= saturateImage(a, satVal)
%% function saturateImage(fnameIn, fnameOut, satVal)
%
%satLim = stretchlim(a, [0, 1-satVal]);
%b=imadjust(a, satLim, [0 1]);
%%b= imadjust(a, [0, 1-satVal], [0 1]);
%
%%-----------------------------------
function b= saturateImage(a, satVal)
% function saturateImage(fnameIn, fnameOut, satVal)
%this assumes 0<a<1

satLim = stretchlim(a(:), [0, 1-satVal]);
for ii = 1:size(a,3)
  a(:,:,ii)=imadjust(a(:,:,ii), satLim, [0 1]);
end
b=a;


%----------------------------------------------------------
function setFigDefaults()


hF =gcf;
hAx = gca;
set(hAx,'TickDir','out')
set(hAx,'Box','off');
set(hF,'Color','w')
set(hAx,'FontSize',14);
%set(hAx,'YTick',[])



%----------------------------------------------------------
function length = getLengthExpt(srInMesh);

nCell = numel(srInMesh);
for ii = 1:nCell
   %LENGTH
   length(ii) = srInMesh{ii}.length;
end


