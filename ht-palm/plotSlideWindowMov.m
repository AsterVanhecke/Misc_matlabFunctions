function plotSlideWindMov(fname, x,y,f,pixSize,varargin)
%function plotSlideWindMov(fname, x,y,f,pixSize,varargin)

sigma=10;
doNorm = false;
gammaVal = 1;
satVal =0.0;
ampScaleFactor = 100;
xLim = [min(x), max(x)];
yLim = [min(y), max(y)];
frameLim = [min(f),max(f)];
tWindow = 500;
deltaT = 50;
srHistVarArg= {};

ii = 1;
while ii <= numel(varargin)
  if strcmp(varargin{ii},'XLim')
    xLim = varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'YLim')
    yLim = varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'FrameLim')
    frameLim = varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Sigma')
    sigma= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Normalise')
     doNorm = true;
     ii = ii +1;
  elseif strcmp(varargin{ii},'Gamma')
    gammaVal= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Saturate')
    satVal= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'AmpScaleFactor')
      ampScaleFactor= varargin{ii+1};
      ii = ii+2;
  elseif strcmp(varargin{ii},'TWindow')
      tWindow= varargin{ii+1};
      ii = ii+2;
  elseif strcmp(varargin{ii},'DT')
      deltaT= varargin{ii+1};
      ii = ii+2;
  elseif strcmp(varargin{ii},'FlipY')
      srHistVarArg={srHistVarArg{:},'FlipY'};
      ii = ii+2;
  else
      ii = ii+1;
  end
end

%filter xy
posOk = (x>=xLim(1) & x<=xLim(2) & y>=yLim(1) & y<=yLim(2));
x = x(posOk);
y = y(posOk);
f = f(posOk);

movieStart = frameLim(1):deltaT:frameLim(2)-tWindow;
movieEnd  = movieStart + tWindow - 1;

for ii = 1:numel(movieStart)
   fMov = (f>=movieStart(ii) & f<=movieEnd(ii));
   xMov = x(fMov);
   yMov = y(fMov);
   
   [srIm] =stormHist2d(xMov,yMov,pixSize,'XLim',xLim,'YLim',yLim,'BlurSigma',sigma,'AmpScaleFactor',ampScaleFactor,srHistVarArg{:});

   if doNorm
      srIm = normSR(srIm,satVal,gammaVal);
   end

   srIm = cast(round(srIm),'uint16');

   if ii == 1
    writeMode = 'overwrite';
   else
    writeMode = 'append';
   end
   imwrite(srIm,fname,'tif','WriteMode',writeMode);
   
end

%--------------------------------------------
function srImOut = normSR(density,satVal,gammaVal)
density = saturateImage(density,satVal);
density = adjustGamma(density,gammaVal);
srImOut = density;



%-----------------------------------------------------------------------------------------------
function imG= adjustGamma(im,gammaVal)

imMax = max(im(:));
%normalise image
imG = ((im/imMax).^gammaVal)*imMax;

%-----------------------------------
function b= saturateImage(a, satVal)
% function saturateImage(fnameIn, fnameOut, satVal)
%outputs assuming uint16
a = (a-min(a(:)))/(max(a(:))-min(a(:)));%normalise
satLim = stretchlim(a(:), [0, 1-satVal]);
a=imadjust(a, satLim, [0 1]);
b=a*65565;

