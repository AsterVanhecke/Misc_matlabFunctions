function plotRingContraction2(srInMeshAll,ringStatName,pixSize,ringRange,tZero,kymographPixSize,lRange,findSpot_xLim,zLim,varargin)

%nCell = numel(srInMeshAll);
%%get the ring only localizations 
%for ii = 1:nCell
%  [l{ii},d{ii},xExt{ii},yExt{ii},z{ii}] = getRingLocs(srInMeshAll,ii,kymographPixSize,lRange,findSpot_xLim,zLim,pixSize,varargin{:});
%  tExpt_min(ii) = srInMeshAll{ii}.tExpt_min;
%end
%
%
%ringCM = mean(ringRange,2);
%poleRingDivide = -550;
%isPolar = [ringCM<poleRingDivide]';
%
%%save them
%save(ringStatName);
load(ringStatName);

diameterLim = 1000;
kk=1;
%for each mid-cell ring 
for ii = 1:nCell
  if ~isPolar(ii) ||isPolar(ii)
    %[zg,ag, bg, alphag] = fitellipseOnAxis2([d{ii},z{ii}]','maxdiameter',diameterLim)
    [zg, rg,initguess] = fitcircle2([d{ii},z{ii}]','maxdiameter',diameterLim,'maxxdisp',75)
    % fit the ring with an ellipse
    %plot an example fit
    hold off;
    plot(d{ii},z{ii},'r.');
    hold all;
    %plotellipse(zg, ag, bg, 0, 'k');
    [xc,yc]=circle(zg(1),zg(2),rg);
    plot(xc,yc,'k');
    [xc,yc]=circle(initguess(1),initguess(2),initguess(3));
    plot(xc,yc,'b-');
    axis equal;
    pause

    fitResult(kk,:) = [kk,isPolar(ii),tExpt_min(ii),rg];
    kk=kk+1
  end
end

%save the radii
save([ringStatName(1:end-4),'ringDiameter.mat'],'fitResult');

%plot the radii vs time


%---------------------------------------
function [x,y,xExt,yExt,z] = getRingLocs(srInMesh,cellNo,kymographPixSize,lRange,findSpot_xLim,zLim,pixSize,varargin)
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

%---------------------------------------
function [xc,yc] = circle(x,y,r)

th = 0:pi/50:2*pi;
xc= r * cos(th) + x;
yc= r * sin(th) + y;

