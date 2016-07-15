function [srIm,xEdge,yEdge] = stormHist2d(x,y,pixSize,varargin)

%vararg defaults
xLim = [min(x) max(x)];
yLim = [min(y) max(y)];
doBlur = false;
blurSigmaPix = 0;
ampScaleFactor = 1;
flipY = false;
doNorm = false;
satVal = 0;
gammaVal = 1;

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
   elseif strcmp(varargin{ii},'BlurSigma')
      doBlur = true;
      blurSigmaRaw = varargin{ii+1};
      blurSigmaPix = blurSigmaRaw/pixSize;
      ii = ii+2;
   elseif strcmp(varargin{ii},'AmpScaleFactor')
      ampScaleFactor= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'FlipY')
      flipY = true;
      ii = ii+1;
   elseif strcmp(varargin{ii},'Gamma')
      gammaVal= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'Normalize')
     doNorm = true;
     satVal = varargin{ii+1};
      ii = ii+2;
   else
      ii = ii+1;
   end
end

[srIm,xEdge,yEdge] = makeSRHist2d(x,y,pixSize,xLim,yLim);
if flipY
  [srIm,yEdge] = doFlipY(srIm,yEdge);
end

srIm = srIm*ampScaleFactor;
if doBlur
   srIm = blurIm(srIm,blurSigmaPix);
end
srIm(srIm<0)=0;
if doNorm ==true
  srIm = saturateImage(srIm,satVal);
end
srIm = adjustGamma(srIm,gammaVal);

%--------------------------------------------------------------
function [srImF,yEdgeF] = doFlipY(srIm,yEdge)

srImF = flipdim(srIm,1);
yEdgeF = yEdge(end:-1:1);
%--------------------------------------------------------------
function [srIm,xEdge,yEdge] = makeSRHist2d(x,y,pixSize,xLim,yLim)

xEdge = xLim(1):pixSize:xLim(2);
yEdge = yLim(1):pixSize:yLim(2);

nx = numel(xEdge);
ny = numel(yEdge);

xInt = interp1(xEdge,1:nx,x,'nearest');
yInt = interp1(yEdge,1:ny,y,'nearest');

%delete any nan points - indicates out of range points
nanPts = isnan(xInt) | isnan(yInt);
xInt(nanPts) = [];
yInt(nanPts) = [];

srIm = accumarray([yInt,xInt],1,[ny,nx]);

%-------------------------------------------------------------
function BW1 = blurIm(BW0,sigmaPsfPix);
hsize=20;

h = fspecial('gaussian',hsize,sigmaPsfPix);
BW1 = imfilter(BW0,h,'replicate');

%-----------------------------------------------------------------------------------------------
function imG= adjustGamma(im,gammaVal)

imMax = max(im(:));
%normalise image
imG = ((im/imMax).^gammaVal)*imMax;

%-----------------------------------
function b= saturateImage(a, satVal)
% function saturateImage(fnameIn, fnameOut, satVal)

a = a/max(a(:));
satLim = stretchlim(a(:), [0, 1-satVal]);
a=imadjust(a, satLim, [0 1]);
b=a;
