function [VOL_blur, hP,isoVal,xEdge, yEdge, zEdge]= Render_3DVol_simple(x,y,z,varargin)

%vararg defaults
pixSize = 20;
sigmaXY = 15;
sigmaZ = 50;
guessIso= true;

xLim = [min(x) max(x)];
yLim = [min(y) max(y)];
zLim = [min(z) max(z)];

if isempty(xLim)|| isempty(yLim)|| isempty(zLim)
   VOL_blur = [];
   hP=[];
   isoVal=0;
   xEdge = [0 1];
   yEdge = [0 1];
   zEdge = [0 1];
else
   %parse varargs
   nn = numel(varargin);
   ii=1;
   while ii<=nn
      if strcmp(varargin{ii},'XLim')
         xLim = varargin{ii+1};
         ii = ii+2;
      elseif strcmp(varargin{ii},'YLim')
         yLim = varargin{ii+1};
         ii = ii+2;
      elseif strcmp(varargin{ii},'ZLim')
         zLim = varargin{ii+1};
         ii = ii+2;
      elseif strcmp(varargin{ii},'PixSize')
         pixSize= varargin{ii+1};
         ii = ii+2;
      elseif strcmp(varargin{ii},'SigmaXY')
         sigmaXY = varargin{ii+1};
         ii = ii+2;
      elseif strcmp(varargin{ii},'SigmaZ')
         sigmaZ = varargin{ii+1};
         ii = ii+2;
       elseif strcmp(varargin{ii},'IsoVal')
         guessIso = false;
         isoVal = varargin{ii+1};
         ii = ii+2;
      else
         ii = ii+1;
      end
   end

   [VOL,xEdge, yEdge, zEdge]= updateVolData(x,y,z,xLim,yLim,zLim,pixSize);

   VOL_blur = updateBlur(VOL,sigmaXY,sigmaZ,pixSize);

   %guess isosurface
   if guessIso
     isoVal = autoIso(VOL_blur);
   end

   hP = plot3d(VOL_blur,isoVal);
end


%---------------------------------------
%helper functions
%-----------------------
function  [VOL,xEdge, yEdge, zEdge]= updateVolData(x,y,z,xLim,yLim,zLim,pixSize)

xEdge = xLim(1):pixSize:xLim(2);
yEdge = yLim(1):pixSize:yLim(2);
zEdge = zLim(1):pixSize:zLim(2);

nx = numel(xEdge);
ny = numel(yEdge);
nz = numel(zEdge);

xInt = interp1(xEdge,1:nx,x,'nearest');
yInt = interp1(yEdge,1:ny,y,'nearest');
zInt = interp1(zEdge,1:nz,z,'nearest');

%delete any nan points - indicates out of range points
nanPts = isnan(xInt) | isnan(yInt) | isnan(zInt);
xInt(nanPts) = [];
yInt(nanPts) = [];
zInt(nanPts) = [];

VOL = accumarray([xInt, yInt, zInt],1,[nx,ny,nz]);

%--------------------------------------
function VOL_blur = updateBlur(VOL,SigmaXY,SigmaZ,pixSize)

blurX = SigmaXY/pixSize;
blurZ = SigmaZ/pixSize;
VOL_blur = gaussf(VOL,[blurX blurX blurZ]);
VOL_blur(VOL_blur<0)=0; %doesn't make sense to have -ive density

%------------------------------------------------
function isoVal = guessIsoVal(VOL_blur)
%guess isosurface
minVol = min(VOL_blur(:));
maxVol = max(VOL_blur(:));
%minVol = max(min(VOL_blur(:)),0);
%maxVol = max(max(VOL_blur(:)),0);
isoVal = (maxVol+minVol)/2;

