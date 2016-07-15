function [srIm,m,n,density,p] = renderData(X,Y,varargin)
%TODO saturate Xpercent before gamma correction
%TODO use a faster gaussian filter
% box = [xstart ystart xwidth yheight]

pixSize=20;
sigma=10;
gammaVal = 1;
satVal =0.0;
zlim = [-600 600];
nz =10;
satVal=0;
isBox = false;
is3D = false;
ii = 1;
while ii <= numel(varargin)
  if strcmp(varargin{ii},'Z')
    is3D = true;
    Z= varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'ZLim')
    zlim = varargin{ii+1};
    ii = ii + 2;
 elseif strcmp(varargin{ii},'Nz')
    nz = varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Pixel Size')
    pixSize = varargin{ii+1};
    ii = ii + 2;
  elseif strcmp(varargin{ii},'Box')
    isBox = true;
    box = varargin{ii+1};
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
  end
end

if is3D ==true
  if isBox ==false
    [srIm ,m,n,density,p] = hist3D(X,Y,Z,pixSize,sigma,gammaVal, zlim,nz,satVal);
  else
    [srIm ,m,n,density,p] = hist3D(X,Y,Z,pixSize,sigma,gammaVal, zlim,nz,satVal,box);
  end
else
  if isBox ==false
    [srIm ,m,n] = hist2D(X,Y,pixSize,sigma,gammaVal,satVal);
  else
    [srIm ,m,n] = hist2D(X,Y,pixSize,sigma,gammaVal,satVal, box);
  end
end

figure;
imagesc(n,m,srIm);
axis equal


%----------------------------------------------------------
function [rgb,m,n,density,p] =hist3D(XPosition,YPosition,ZPosition,pixSize,sigma,gammaVal,zlim,nz, satVal,box)

minC = 0;
maxC =1;

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

if ~exist('zlim','var')
  minZ = min(ZPosition);
  maxZ = max(ZPosition);
else
  minZ = zlim(1);
  maxZ = zlim(2);
end

%remove out of bounds data
isInBounds = XPosition > minX & XPosition < maxX ...
              & YPosition > minY & YPosition < maxY ...
              & ZPosition > minZ & ZPosition < maxZ;

XPosition = XPosition(isInBounds);
YPosition = YPosition(isInBounds);
ZPosition = ZPosition(isInBounds);

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
%ZPosition=evalin('base',handles.varz);
% [minZ maxZ]=getZbounds(handles);
notrials = 1; % SET NOTRIALS TO 20 FOR JITTERING
p=minZ:(maxX-minZ)/nz:maxZ;
NN = numel(XPosition);
%amount = (evalin('base',handles.sigmax)+evalin('base',handles.sigmay))/2;
amount = (sx+sy)/2;
amount = amount*0; % TO AVOID JITTERING
%density = zeros(res,res,nz);
density = zeros(ny,nx,nz);
%setProgress(handles,0);
%pxx=(maxX-minX)/res; pxy=(maxY-minY)/res; pxz=(maxZ-minZ)/res;
pxx=pixSize; pxy=pixSize; pxz=(maxZ-minZ)/nz;
pxVol=pxx * pxy * pxz;
for i=1:notrials
%    RR = [round((res-1)*[...
%        (XPosition(subset)+randn(NN,1).*amount(subset)-minX)/(maxX-minX) ...
%        (YPosition(subset)+randn(NN,1).*amount(subset)-minY)/(maxY-minY) ])+1 ...
%          round(nz*(ZPosition(subset)+randn(NN,1).*amount(subset)-minZ)/(maxZ-minZ))];
%    RRok = all(RR(:,1:2)<=res,2) & all(RR>=1,2) & RR(:,3)<=nz;
    RR = [...
          round((ny-1)*(YPosition(subset)+randn(NN,1).*amount(subset)-minY)/(maxY-minY))+1 ...
          round((nx-1)*(XPosition(subset)+randn(NN,1).*amount(subset)-minX)/(maxX-minX))+1, ...
          round((nz-1)*(ZPosition(subset)+randn(NN,1).*amount(subset)-minZ)/(maxZ-minZ))+1];
    RRok = logical(ones(size(RR,1),1));
    %d = accumarray(RR(RRok,:),1,[res,res,nz])/pxVol;
    d = accumarray(RR(RRok,:),1,[ny,nx,nz])/pxVol;
    if i==1
        density = d;
    else
        density = density + d;
    end

    %X = repmat(n,res,1); Y = repmat(m',1,res);
    density = density/notrials;
    
%    setProgress(handles,1/notrials)
%    drawnow
end


%%TEST EXPORT 3D STACK
%%<DEBUG
%maxD = max(density(:));
%for ii = 1:nz
%  if ii == 1
%    writeMode = 'overwrite';
%  else
%    writeMode = 'append';
%  end
%  dIm = im2uint8(density(:,:,ii)/maxD);
%  imwrite(dIm,'test.tif','tif','WriteMode',writeMode);
%end
%%DEBUG>

%sigma = str2double(get(handles.sigma,'String'));

%TODO use a faster gauss filter here
sPix = sigma/pxx;
gWindow = ceil(5*sPix);
gKern = fspecial('gaussian',gWindow, sPix);
dMax = max(density(:));
density = density/max(density(:));
for ii = 1:size(density,3)
  density(:,:,ii) = imfilter(density(:,:,ii),gKern,'replicate');
end
% ADJUST GAMMA HERE
density(density<0)=0;
density = saturateImage(density,satVal);
density = adjustGamma(density,gammaVal);

%convert density to dipimage
density = dip_image(density);
HSV = joinchannels('HSV',zz(density,'corner')/(nz+1)*2*pi,1+newim(size(density)),density);
RGB = colorspace(HSV,'RGB');

a = .5; b=1; sbar = 1;

rgb = zeros(size(density,2),size(density,1),3);
for i=1:size(density,3)
    for c=1:3
        rgb(:,:,c) = rgb(:,:,c).*double(1-a*squeeze(density(:,:,i-1))) + double(b*squeeze(RGB{c}(:,:,i-1)*density(:,:,i-1)));
    end
end

adjustBrightness= @(im,minC,maxC)  max(0,min((im/max(im(:)) - minC) / (maxC-minC),1));
rgb = adjustBrightness(rgb,minC,maxC);


%rgb(rgb<0)=0; 
%rgb(rgb>1)=1;

%sz = [res/4 res/32]; 
%szt = [res/4 3]; 
%cbx = 10; cby = res-sz(2)-10*(res/256);
%szs = [round(sbar/pxx) res/128];
%sbx = res-5*(res/256)-szs(1); sby = res-szs(2)-20*(res/256);
%
%maxc = max(rgb(:));
%
%tickss = colorspace(...
%    joinchannels('HSV',...
%    xx(szt,'corner')/szt(1)*2*pi,...
%    newim(szt)+1,...
%    newim(szt)+maxc),'RGB');
%for i=0:10:(szt(1)-1)
%    tickss{:}(i,:)=1;
%end
%
%colb = colorspace(...
%    joinchannels('HSV',...
%    xx(sz,'corner')/sz(1)*2*pi,...
%    newim(sz)+.5,...
%    (1-yy(sz,'corner')/sz(2))*maxc),'RGB');
%
%rgb(cby+(1:sz(2)),cbx+(1:sz(1)),:)=cat(3,double(colb{1}),double(colb{2}),double(colb{3}));
%rgb(cby+(1:szt(2))-3,cbx+(1:szt(1)),:)=cat(3,double(tickss{1}),double(tickss{2}),double(tickss{3}));
%rgb(sbx+(1:szs(2)),sby+(1:szs(1)),:)=maxc;


%text(minX+cbx*pxx,minY+pxy*(cby-5*(res/256)),num2str(minZ),'Color','w', 'HorizontalAlignment','left')
%text(minX+(cbx+sz(1))*pxx,minY+pxy*(cby-5*(res/256)),num2str(maxZ),'Color','w', 'HorizontalAlignment','right')
%text(minX+pxx*(sby+szs(1)/2),minY+pxy*(sbx-5*(res/256)),[num2str(100) ' nm'],'Color','w', 'HorizontalAlignment','center')
%text(minY+pxy*(cbx+sz(1)),minX+pxx*(cby-5),num2str(maxZ),'Color','w', 'HorizontalAlignment','right')
        
%-----------------------------------------------------------------------------------------------
function [density,m,n] =hist2D(XPosition,YPosition,pixSize,sigma,gammaVal,satVal,box)

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

%TODO use a faster gauss filter here
sPix = sigma/pxx;
gWindow = ceil(5*sPix);
gKern = fspecial('gaussian',gWindow, sPix);
dMax = max(density(:));
density = density/max(density(:));
density = imfilter(density,gKern,'replicate');

%density = gaussf(density/max(density(:)),[sigma/pxx sigma/pxy 0]);
% ADJUST GAMMA HERE
density(density<0)=0;
density = saturateImage(density,satVal);
density = adjustGamma(density,gammaVal);

%-----------------------------------------------------------------------------------------------
function imG= adjustGamma(im,gammaVal)

imMax = max(im(:));
%normalise image
imG = ((im/imMax).^gammaVal)*imMax;

%-----------------------------------
function b= saturateImage(a, satVal)
% function saturateImage(fnameIn, fnameOut, satVal)
%this assumes 0<a<1

satLim = stretchlim(a(:), [0, 1-satVal]);
for ii = 1:size(a,3)
  a(:,:,ii)=imadjust(a(:,:,ii), satLim, [0 1]);
end
b=a;
