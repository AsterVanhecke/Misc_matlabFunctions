function [stdR_nm ,tExpt_min,rCM, rcmBW, sR_BW , goodIm]= calculateRingRadius(srInMesh,movieName,ringStatName,pixSize)

minPix = 5;

info = imfinfo(movieName);
nFrame = numel(info);
for ii = 1:nFrame
   im0(:,:,ii) = imread(movieName,'Index',ii,'Info',info);
end

%normalise im0
norm = @(x) (x - min(x(:)))./(max(x(:)) - min(x(:)));

im0 = double(im0);
im0 = norm(im0);

% use otsus method to segment the images
lvl = graythresh(im0);
imBW_All = (im0>lvl);

%discard frames that have fewer than minPix counts 
for ii = 1:nFrame
   imCur = imBW_All(:,:,ii);
   goodIm(ii) = sum(imCur(:)) > minPix;
end

%max individually segmented image too
imBW = imBW_All;%initialise memory
for ii = 1:nFrame
   imC = im0(:,:,ii);
   imBW(:,:,ii) = imC>graythresh(imC);%assign
end

%hold off;
%for ii = 1:nFrame
%   subplot(1,3,1);
%   imshow(imBW(:,:,ii));
%   subplot(1,3,2);
%   imagesc(im0(:,:,ii));
%   axis equal;
%   subplot(1,3,3);
%   imC = im0(:,:,ii);
%   imshow(imC>graythresh(imC));
%   if goodIm(ii)
%      title('Good');
%   else 
%      title('bad');
%   end
%   pause
%end

mean_nm = zeros(nFrame,2);
std_nm = zeros(nFrame,2);
tExpt_min = zeros(nFrame,1);
maxRad = floor(min(info(1).Width,info(1).Height)/2);
intRadial = zeros(nFrame,maxRad);
for ii = 1:nFrame
   im = im0(:,:,ii);
   
   [ux,uy,sx,sy] = getImMoment(im);

   [rCM_pix] = getRadialProfile(im,ux,uy,maxRad);
   rCM(ii) = rCM_pix*pixSize;

   mean_nm(ii,:) = [ux,uy].*pixSize;
   std_nm(ii,:) = [sx,sy].*pixSize;
   tExpt_min(ii) = srInMesh{ii}.tExpt_min;
   
   [uxBW,uyBW,sxBW,syBW] = getImMoment(imBW(:,:,ii));
   [rcmBW_pix] = getRadialProfile(imBW(:,:,ii),uxBW,uyBW,maxRad);
   rcmBW(ii) = rcmBW_pix*pixSize;
   meanBW_nm(ii,:) = [uxBW,uyBW].*pixSize;
   stdBW_nm(ii,:) = [sxBW,syBW].*pixSize;
   sR_BW(ii) = sqrt(sxBW.^2+syBW.^2).*pixSize;



   %ii
   %sR = sqrt(sx.^2+sy.^2);
   %[xc1,yc1] = circle(ux,uy,sR);
   %[xc2,yc2] = circle(ux,uy,rCM_pix);
   %[xc3,yc3] = circle(uxBW,uyBW,sR_BW);
   %[xc4,yc4] = circle(uxBW,uyBW,rcmBW_pix);
   %hold off;
   %imagesc(im);
   %colormap('gray');
   %hold all;
   %plot(ux,uy,'bx')
   %plot(xc1,yc1,'y-')
   %plot(xc2,yc2,'r-')
   %plot(xc3,yc3,'b-')
   %plot(xc4,yc4,'g-')
   %legend('Centre','StDev','RadialCM','Stdev BW', 'RadCM BW');
   %if goodIm(ii)
   %   title('Good');
   %else 
   %   title('bad');
   %end
   %pause
   
end

stdR_nm = sqrt(sum(std_nm.^2,2));
%
%h1 = figure;
%hold all;
%plot(tExpt_min,stdR_nm,'k.')
%plot(tExpt_min,rCM,'b.')
%xlabel('Time, min');
%ylabel('Radial Std Dev, nm');
%legend('Stdev','R CM','R Max');
%saveas(h1,[ringStatName,'_stdevPlot1.fig']);

%---------------------------------------
function [xc,yc] = circle(x,y,r)

th = 0:pi/50:2*pi;
xc= r * cos(th) + x;
yc= r * sin(th) + y;

%---------------------------------------
function [ux,uy,sx,sy] = getImMoment(im);

iX = sum(im,1);
iY = sum(im,2);
x = 1:size(im,2);
y = [1:size(im,1)]';
iXT = sum(iX);
iYT = sum(iY);
ux = sum(x.*iX)/iXT;
uy = sum(y.*iY)/iYT;
sx = sqrt(sum(iX.*(x-ux).^2)/iXT);
sy = sqrt(sum(iY.*(y-uy).^2)/iYT);
