function autoAlignChannel(acqFolder,acqConfigFile, saveFolder, tformSaveName,showFits)

TFORM_TYPE = 'lwm';
IM_EDGE=5;
FIT_WINDOW=30;

if ~exist('showFits','var')
   showFits=false;
end


% import the xml file
config = xml_read([acqFolder,filesep(),acqConfigFile]);% using xml_read from xml_io_tools_2010_11_0
acqFileBaseName = config.fileBaseName_;
acqMetadataFname = [acqFolder,filesep(),acqFileBaseName,'_htpalmMetadata.xml'];
acqMetadata = xml_read(acqMetadataFname);

folderNoStub = '_1';
imStub = '_MMImages.ome.tif';
imMetadataStub = '_MMImages_metadata.txt';

% get the list of fl images
% and the list of phPre images
nAcq = acqMetadata.nAcq_;
for ii = 1:nAcq
   flName{ii} = [acqMetadata.acqMetadataList_.acqMetadata(ii).acqNameFl(1).CONTENT];
   flImName{ii} = [acqFolder,filesep(),flName{ii},folderNoStub,filesep(),flName{ii},imStub];
   flImMeta{ii} = [acqFolder,filesep(),flName{ii},folderNoStub,flName{ii},imMetadataStub];

   phName{ii} = [acqMetadata.acqMetadataList_.acqMetadata(ii).acqNamePhPost_];
   phImName{ii} = [acqFolder,filesep(),phName{ii},folderNoStub,filesep(),phName{ii},imStub];
   phImMeta{ii} = [acqFolder,filesep(),phName{ii},folderNoStub,phName{ii},imMetadataStub];
end

% get the pixel size for each camera
imSummaryPh = readImSummaryMM(phImName{1});
imSummaryFl = readImSummaryMM(flImName{1});
pixUmPh = imSummaryPh.PixelSize_um;
pixUmFl = imSummaryFl.PixelSize_um;
roiPixPh = imSummaryPh.ROI;
roiPixFl = imSummaryFl.ROI;

hF = figure;
% get the initial centre and offsets for each channel using non-cpselect based tool
for ii = 1:3
   %get the ph im bead centre
   phIm = imread(phImName{ii});
   phIm = imadjust(phIm, stretchlim(phIm,0),[0 1]);
   set(0,'CurrentFigure',hF);
   pos = getBeadPosManual(phIm,FIT_WINDOW);
   phPosPix(ii,:) = pos;

   %get the fl im bead centre
   flIm = imread(flImName{ii});
   flIm = imadjust(flIm, stretchlim(flIm,0),[0 1]);%invert the phase contrast image
   set(0,'CurrentFigure',hF);
   pos = getBeadPosManual(flIm,FIT_WINDOW);
   flPosPix(ii,:) = pos;
end
close(hF);

% and calculate the list of offsets from initial position
stepXUm = config.mosaicStepSizeX_;
stepYUm = config.mosaicStepSizeY_;

initPosListPhPix  = getSpiralPos(phPosPix,stepXUm/pixUmPh,stepYUm/pixUmPh,nAcq);
initPosListFlPix  = getSpiralPos(flPosPix,stepXUm/pixUmFl,stepYUm/pixUmFl,nAcq);

%trim the list to only positions within the FOV
insideFovPh = ( initPosListPhPix(:,1)> IM_EDGE &  initPosListPhPix(:,1) < (roiPixPh(3)-IM_EDGE)...
                  & initPosListPhPix(:,2)> IM_EDGE &  initPosListPhPix(:,2) < (roiPixPh(4)-IM_EDGE));
insideFovFl = ( initPosListFlPix(:,1)> IM_EDGE &  initPosListFlPix(:,1) < (roiPixFl(3)-IM_EDGE)...
                  & initPosListFlPix(:,2)> IM_EDGE & initPosListFlPix(:,2) < (roiPixFl(4)-IM_EDGE));

%for each bead inside the FOV, fit the position of the bead
if showFits==true
   hF = figure;
end

jj = 1;
for ii = 1:nAcq
   if insideFovPh(ii) && insideFovFl(ii)
      %get the ph im bead centre
      phIm = imread(phImName{ii});
      phIm = imadjust(phIm, stretchlim(phIm,0),[0 1]);
      x = initPosListPhPix(ii,1);
      y = initPosListPhPix(ii,2);
      pos = getPosGauss2d_PH(phIm,[x,y],FIT_WINDOW);
      %pos = getPeakPosCentroidPH(phIm,[x,y],FIT_WINDOW);
      posPixPh(jj,:) = pos;

      if showFits==true
         set(0,'CurrentFigure',hF);
         hold off;
         imshow(phIm);
         hold all;
         plot(pos(1),pos(2),'x');
         pause
      end

      %get the fl im bead centre
      flIm = imread(flImName{ii});
      flIm = imadjust(flIm, stretchlim(flIm,0),[0 1]);%invert the phase contrast image
      x = initPosListFlPix(ii,1);
      y = initPosListFlPix(ii,2);
       [pos amplitude sigma_x sigma_y]= getPosGauss2d(flIm,[x,y],FIT_WINDOW);
      posPixFl(jj,:) = pos;
      sigmaPixFl(jj,:) = [sigma_x,sigma_y];

      if showFits==true
         set(0,'CurrentFigure',hF);
         hold off;
         imshow(flIm);
         hold all;
         plot(pos(1),pos(2),'x');
         pause
      end

      jj=jj+1;
   end
end

%% correct for xy wobble
%%posPixFl = correctCalWobble(posPixFl,sigmaPixFl,wobbleFile,srPixSize);

%apply the roi offset to each coord set
posPixFl_ROI(:,1) = posPixFl(:,1) + roiPixFl(1);
posPixFl_ROI(:,2) = posPixFl(:,2) + roiPixFl(2);
posPixPh_ROI(:,1) = posPixPh(:,1) + roiPixPh(1);
posPixPh_ROI(:,2) = posPixPh(:,2) + roiPixPh(2);

%use cp2tform to generate TFORM
dblTform_phToFl = getDblTform(posPixPh_ROI,posPixFl_ROI,TFORM_TYPE);
save(tformSaveName,'dblTform_phToFl','posPixFl_ROI','posPixPh_ROI','posPixFl', 'posPixPh');

[XIn] = dblTformInv(dblTform_phToFl,posPixFl_ROI);
hold off;
plot(posPixPh_ROI(:,1),posPixPh_ROI(:,2),'kx');
hold all;
plot(posPixFl_ROI(:,1),posPixFl_ROI(:,2),'bx');
plot(XIn(:,1),XIn(:,2),'rx');
%--------------------------------------------
function pos = getBeadPosManual(im,windowSize)

hold off;
imshow(im);
[x,y] = ginput(1);
hold all;
plot(x,y,'gx');
pos = getPosGauss2d(im,[x,y],windowSize);
plot(pos(1),pos(2),'rx');

%--------------------------------------------
function posList= getSpiralPos(first3Pos,stepX,stepY,nAcq);

[dX X0 theta] = getSpiralParam(first3Pos);
if theta>0
   R = [0,-1;
        1 0];
else
   R = -1*[0,-1;
           1 0];
end

%calculate all the coordinates
posList = zeros(nAcq,2);
posList(1,:) = X0;

%integer coordinates to work out when to rotate
dXint_=1; 
dYint_=0; 
xInt_=0;
yInt_=0;

for ii = 2:nAcq
   xInt_ = xInt_ + dXint_ ;
   yInt_ = yInt_ + dYint_;
   posList(ii,:) = posList(ii-1,:)+dX;

   %update dy dx
   if ( (xInt_==yInt_) || (xInt_ < 0 && xInt_ == -yInt_) || (xInt_ > 0 && xInt_ == 1-yInt_) ) %if at edge
      % apply a 90deg rotation to (dXint, dYint)
      %(dX') = (0 -1)(dX)
      %(dY')   (1  0)(dY)
      t = dXint_;
      dXint_ = -dYint_;
      dYint_ = t;

      %apply a theta (either +,- 90) rotation to dX
      dX = (R*dX')';
   end
end
%
%offsetList = zeros(nFov,2);
%dXint_=1; 
%dYint_=0; 
%xInt_=0;
%yInt_=0; 
%nCur_=0;
%for ii = 2:nFov
%   [x,y, dXint_,dYint_,xInt_,yInt_] = gotoNextFov( dXint_,dYint_,xInt_,yInt_,stepX,stepY);
%   offsetList(ii,:) = [x,y];
%end

%--------------------------------------------
function [dX X0 theta] = getSpiralParam(first3Pos)
X0 = first3Pos(1,:);
dX = first3Pos(2,:)-first3Pos(1,:);
dX2 = first3Pos(3,:)-first3Pos(2,:);
% find angle of dX - gives tilt of system
% find angle between dX1, dX2 - gives direction of rotation of system
thetaTmp = atan2(dX2(2),dX2(1)) -  atan2(dX(2),dX(1));
if thetaTmp>0
   theta = pi/2;
else
   theta = -pi/2;
end


%--------------------------------------------
function [x_,y_, dXint_,dYint_,xInt_,yInt_] = gotoNextFov( dXint_,dYint_,xInt_,yInt_,xStep_,yStep_)

xStart_ = 0;
yStart_ = 0;

xInt_ = xInt_ + dXint_ ;
yInt_ = yInt_ + dYint_;
x_ = xInt_*xStep_+xStart_;
y_ = yInt_*yStep_+yStart_;

%update dy dx
if ( (xInt_==yInt_) || (xInt_ < 0 && xInt_ == -yInt_) || (xInt_ > 0 && xInt_ == 1-yInt_) ) %if at edge
      % apply a 90deg rotation to (dX, dY)
      %(dX') = (0 -1)(dX)
      %(dY')   (1  0)(dY)
      t = dXint_;
      dXint_ = -dYint_;
      dYint_ = t;
end
  
%-------------------------------------------
function posPixFlOut = correctCalWobble(posPixFl,sigmaPixFl,wobbleFile,srPixSize);
% correct for xy wobble

calData = importdata(wobbleFile);

%convert sigma to nm
sigmaXY_nm = [sigmaPixFl(:,1)*srPixSize(1), sigmaPixFl(:,2)*srPixSize(2)];
zPos = getZPosSimple(sigmaXY_nm, calData(:,1:3));
avgZ = median(zPos);

xPosNm = [posPixFl(:,1)*srPixSize(1)];
yPosNm = [posPixFl(:,2)*srPixSize(2)];
zPosNm = ones(size(xPosNm))*avgZ;

zScaleFactor =1 ;%ignore z rescaling here since we will throw away z!
[xC,yC,zC] = correct3D(xPosNm,yPosNm,zPosNm,wobbleFile,zScaleFactor);

posPixFlOut = [xC(:)/srPixSize(1), yC(:)/srPixSize(2)];
