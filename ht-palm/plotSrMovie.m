function plotSrMovie(srInMesh,ringRange,gammaVal,satVal,movieFolder, movieName, bactList,zLim,pixSize,tWindow,deltaT,frameLim,sigma,varargin)

if ~exist(movieFolder,'dir')
   mkdir(movieFolder);
end

xMin = -4000;xMax = 4000;
yMin = -400;yMax = 400;
zMin = zLim(1);zMax = zLim(2);
ampScaleFactor = 100;
doSortCells=false;
doFlipZ=false;
zMovVarArg={};
ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'SortCells')
    doSortCells=true;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'FlipZ')
    doFlipZ=true;
    zMovVarArg = {zMovVarArg{:},'FlipY'};
    ii = ii + 1;
  else
    ii = ii + 1;
  end
end

m=xMin:pixSize:xMax;
n=yMin:pixSize:yMax;
p = zMin:pixSize:zMax;


cellLength = getCellLength(srInMesh);
nCell = numel(cellLength);

nFullCell = 0;
if doSortCells
   [dum, order]=sort(cellLength(:,1));
else
   order = 1:nCell;
end

for kk = bactList
  ii = order(kk);
  
  srData = srInMesh{ii}.localizations;
  mesh   = srInMesh{ii}.meshAligned;

  %only show points inside ymin rand and x slice lRange
  x=srInMesh{ii}.l;
  y=srInMesh{ii}.d;
  zcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','position in sample space in z dimension');
  z=srInMesh{ii}.localizations.data(:,zcol);
  fcol = findSRField(srInMesh{ii}.localizations.parInfo,'semantic','frame number');
  f=srInMesh{ii}.localizations.data(:,fcol);

  inCell = (y>yMin & y<yMax);
  x=x(inCell);
  y=y(inCell);
  z=z(inCell);
  f=f(inCell);

  %recentre and flip the data along CofM
  [x,y,mesh] = alignCM(x,y,mesh);

   
  %exclude points outside of the lRange
  inCell = (x>ringRange(kk,1) & x<ringRange(kk,2));
  xCross=x(inCell);
  yCross=y(inCell);
  zCross=z(inCell);
  fCross=f(inCell);

  %plot xy movie
  fname1 = fullfile(movieFolder,[movieName,'_ld_',num2str(kk),'.tif']);
  plotSlideWindowMov(fname1, x,y,f,pixSize,'XLim',[xMin,xMax],'YLim',[yMin,yMax],'FrameLim',frameLim,...
               'Sigma',sigma,'AmpScaleFactor',ampScaleFactor,'TWindow',tWindow,'DT',deltaT);

  %plot yz movie 
  fname1 = fullfile(movieFolder,[movieName,'_dz_',num2str(kk),'.tif']);
  plotSlideWindowMov(fname1, yCross,zCross,fCross,pixSize,'XLim',[yMin,yMax],'YLim',[zMin,zMax],'FrameLim',frameLim,...
               'Sigma',sigma,'AmpScaleFactor',ampScaleFactor,'TWindow',tWindow,'DT',deltaT,zMovVarArg{:});
  %[srIm] =stormHist2d(y,z,pixSize,'XLim',[yMin,yMax],'YLim',[zMin,zMax],'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor);

  %save the localizations
  fname1 = fullfile(movieFolder,[movieName,'_dz_',num2str(kk),'.ringLocs.mat']);
  save(fname1,'xCross','yCross','zCross','fCross');
end
