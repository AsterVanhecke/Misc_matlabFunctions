function [xDriftNm,yDriftNm] =phFlCombine3(srFile,meshFile, srInMeshName, tformFname, srPixSize, phPixSize,roiPixFl, roiPixPh,units,phPathPre,phPathPost,zScaleFactor,nFrameFl,tExpt_min,xyWobbleFile)
try 
   m = mTrack();

   %Work out the right file names

   imParamPh.pixSize = phPixSize;
   imParamPh.ROI     = roiPixPh;
   %account for zeroindex MM
   imParamPh.ROI(1:2) = imParamPh.ROI(1:2) +1;

   imParamFl.pixSize = srPixSize;
   imParamFl.ROI     = roiPixFl;
   %account for zeroindex MM
   imParamFl.ROI(1:2) = imParamFl.ROI(1:2) +1;
   %<DEBUG
   %imParamFl.ROI(1:2) = [1 1]
   %DEBUG>

   % load the SR analysis
   srData = readRapidStorm(srFile);
   xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
   ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
   zcol = findSRField(srData.parInfo,'semantic','position in sample space in Z');
   fcol = findSRField(srData.parInfo,'semantic','frame number');
   photCol =  findSRField(srData.parInfo,'semantic','emission strength');
   if zcol < 0
     srIs3D =false;
   else
     srIs3D = true;
   end

   srDataOld = srData;
   %apply 3d correction factors
   if srIs3D
      srData = apply3DCorrect(srData,xyWobbleFile,zScaleFactor);
   end

   %transform to PH coordinate system
   [srData xy]= transformToPH(srData,srPixSize,phPixSize,tformFname,imParamPh,imParamFl);

   %%%<TEMP HACK TO FIX BROKEN TFORM
   %%generateTformFromData(xy,phPathPre,srPixSize,phPixSize,tformFname,imParamPh,imParamFl);
   %xDriftNm = NaN;
   %yDriftNm = NaN;
   %xyRoiAll ={};
   %roiDataFname = '150206roiData.mat';
   %if exist(roiDataFname,'file')
   %  load(roiDataFname,'xyRoiAll');
   %end
   %xyRoiAll ={xyRoiAll{:},xy};
   %save(roiDataFname,'xyRoiAll');
   %return 
   %%%TEMP HACK TO FIX BROKEN TFORM>

   %linear drift correction based on drift in phase contrast images
   srDataO = srData;
   [srData, xDriftNm,yDriftNm] = driftCorrectPH(srData,phPathPre,phPathPost,nFrameFl,phPixSize);

   %<DEBUG
   %hold off;
   %plot(srDataO.data(:,xcol),srDataO.data(:,ycol),'r.','MarkerSize',6);
   %hold all;
   %plot(srData.data(:,xcol),srData.data(:,ycol),'k.','MarkerSize',6);
   %axis equal;
   %pause
   %DEBUG>

   %load the cell list
   meshData = load(meshFile);
   cellList = meshData.cellList;

   %%<DEBUG
   % TEMP: plot the overlays
   %phIm = imread(phPathPre);
   %phIm = imadjust(phIm, stretchlim(phIm,0),[0 1]);%invert the phase contrast image
   %hold off;
   %imshow(phIm);
   %hold all;
   %dispcellall(cellList,phIm);
   %hold all;
   %plot(xy(:,1),xy(:,2),'rx','MarkerSize',2);
   %pause
   %%DEBUG>

   cellList= changeUnitSystem(cellList, phPixSize(1),units);

   %% delete cells whos box intersects image edge or is outside image
   datatmp = load(tformFname,'dblTform_phToFl');
   dblTform_phToFl= datatmp.dblTform_phToFl;
   cellList = deleteOutsideCell(cellList,imParamFl.ROI,dblTform_phToFl);

   %%<DEBUG
   %%plot the image and localisations
   %
   %gammaVal = 1;
   %satVal = 0.001;
   %%satVal = 0.0001;
   %%satVal = 0.0;
   %pixSize=20;sigma=10;
   %xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
   %ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
   %zcol = findSRField(srData.parInfo,'semantic','position in sample space in Z');
   %x=srData.data(:,xcol);y=srData.data(:,ycol);z=srData.data(:,zcol);
   %
   %%[srIm,m,n] = renderData(x,y,'Z',z,'Pixel Size',pixSize,'Sigma',sigma,'Gamma',gammaVal,'Saturate',satVal);
   %[srIm,m,n] = renderData(x,y,'Pixel Size',pixSize,'Sigma',sigma,'Gamma',gammaVal,'Saturate',satVal);
   %hold off;
   %colormap('hot');
   %imagesc(n,m,srIm);axis equal
   %hold all;
   %%plot(xyDataPH(:,1),xyDataPH(:,2),'rx');
   %dispcellallSR(cellList,phIm);
   %pause
   %%DEBUG>

   % put SR points in PH meshes (use spotfinder routines)
   cellList = associatePhSr(cellList, srData);

   cellList = addExptTime(cellList,tExpt_min);

catch ME
   %if anything goes wrong skip this cell, assign empty dataset
   cellList = cell(1);
   cellList{1} = cell(1);
   xDriftNm = NaN;
   yDriftNm = NaN;
   fprintf('\n**************************************************\n');
   fprintf(  '******* File skipped for following reason: *******\n');
   fprintf('%s\n',getReport(ME));
   fprintf(  '**************************************************\n');
end
save(srInMeshName,'cellList','srData');

%%% NB HOW to plot the meshes:
%%%mesh = cellList{frame}{cell}.mesh;
%%%plot(mesh(:,1),mesh(:,2))%left side of cell
%%%plot(mesh(:,3),mesh(:,4));%right side of cell
%%%plot(mesh(:,[1 3])',mesh(:,[2 4])')%segments
%----------------------------------------
function doSegmentation(phFolder,saveName,algFname);

m=mTrack;

m.loadimages(1,phFolder);
m.loadparam(algFname);

range=[1 1];
mode=3;
lst=[];
addsig=[];
addas={};
savefile=saveName;
fsave=[];
saveselect=0;
processregion=[];

m.process(range,mode,lst,addsig,addas,savefile,fsave,saveselect,processregion);

%---------------------------------------------------------------------
function cellListOut = changeUnitSystem(cellList, pixSize,units)
%pixSize pix
%units: 'pixels', 'nanometres'

nFrame = numel(cellList);
%for each frame
for ii = 1:nFrame
  nCell = numel(cellList{ii});
  %for each cell
  for jj = 1:nCell
    if ~isempty(cellList{ii}{jj}) 
      if ~isfield(cellList{ii}{jj},'length')
         %correct for weird MT bug where a box gets assigned to a cell,
         % but theres no cell inside - delete these
         cellList{ii}{jj}=[]; 
      else
         %add the coordSystemParam structure
         cellList{ii}{jj}.pixSize  = pixSize;
         cellList{ii}{jj}.units    = units;
         %transform the x,y coordinates:(offset is applied to these)
         % box,mesh
         cellList{ii}{jj}.box = cellList{ii}{jj}.box*pixSize;;
         cellList{ii}{jj}.mesh = cellList{ii}{jj}.mesh*pixSize;
         % length/area/volume
         % steplength/steparea/stepvolume 
         cellList{ii}{jj}.length	= cellList{ii}{jj}.length*pixSize;
         cellList{ii}{jj}.area	=cellList{ii}{jj}.area*pixSize;
         cellList{ii}{jj}.volume	=cellList{ii}{jj}.volume*pixSize;
         cellList{ii}{jj}.steplength	= cellList{ii}{jj}.steplength*pixSize;
         cellList{ii}{jj}.steparea	= cellList{ii}{jj}.steparea	*pixSize;
         cellList{ii}{jj}.stepvolume	=cellList{ii}{jj}.stepvolume*pixSize;
      end
    end
  end
end

cellListOut = cellList;

%----------------------------------------------------------------
function dispcellallSR(cellList,phIm)

%imshow(phIm);
hold all;

nFrame = numel(cellList);
%for each frame
for ii = 1:nFrame
  nCell = numel(cellList{ii});
  %for each cell
  for jj = 1:nCell
    if ~isempty(cellList{ii}{jj})
      pixSize = cellList{ii}{jj}.pixSize;
      mesh = cellList{ii}{jj}.mesh;
      %mesh = mesh/pixSize;
      plot(mesh(:,1),mesh(:,2),'g-')%left side of cell
      plot(mesh(:,3),mesh(:,4),'g-');%right side of cell
      %plot(mesh(:,[1 3])',mesh(:,[2 4])','k-')%segments
    end
  end
end
%-----------------------------------------------------
function cellListOut = deleteOutsideCell(cellList,RoiFlBox,dblTform_phToFl)
% delete cell whose box intersects image edge or are outside image

%convert to pixel edges
%make Fl pixel ROI box
RoiFl = box2Vertex(RoiFlBox);

%transform to Ph pixel
m = mTrack();
RoiFlPH = m.dblTformInv(dblTform_phToFl,RoiFl);

%transform to Ph nano
pixSize= cellList{1}{1}.pixSize;
RoiFlPH = RoiFlPH*pixSize;

%%<DEBUG
%hold off
%plot([RoiFl(:,1);RoiFl(1,1)],[RoiFl(:,2);RoiFl(1,2)],'r-');
%hold all
%plot([RoiFlPH(:,1);RoiFlPH(1,1)],[RoiFlPH(:,2);RoiFlPH(1,2)],'k-');
%%DEBUG>

% For each cell, check if box is entirely inside ROI, else delete

nFrame = numel(cellList);
%for each frame
for ii = 1:nFrame
  nCell = numel(cellList{ii});
  %for each cell
  for jj = 1:nCell
    if ~isempty(cellList{ii}{jj})
      box= cellList{ii}{jj}.box;
      boxV = box2Vertex(box);
      boxIN = inpolygon(boxV(:,1),boxV(:,2),RoiFlPH(:,1),RoiFlPH(:,2));
      if ~all(boxIN == 1)
        cellList{ii}{jj}=[];
      end
    end
  end
end
cellListOut = cellList;

%------------------------------------------
function boxV = box2Vertex(box)

boxV =[box(1)         , box(2);
        box(1)+box(3)-1, box(2);
        box(1)+box(3)-1, box(2)+box(4)-1;
        box(1)         , box(2)+box(4)-1];
%-------------------------------------------------------
function xyDataPH = nanOutsidePh(xyDataPH,RoiPhBox);
%change points at image edges which may not have been sucessfully transformed to nan
xMax = RoiPhBox(3);
yMax = RoiPhBox(4);

badPts = find(xyDataPH(:,1) >xMax | xyDataPH(:,2) >yMax);
xyDataPH(badPts,:) = NaN;
%-------------------------------------------------------
function srDataOut = apply3DCorrect(srData,wobbleFile,zScaleFactor)
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
zcol = findSRField(srData.parInfo,'semantic','position in sample space in Z');

x = srData.data(:,xcol);
y = srData.data(:,ycol);
z = srData.data(:,zcol);
[xC,yC,zC] = correct3D(x,y,z,wobbleFile,zScaleFactor);

srData.data(:,xcol) = xC; 
srData.data(:,ycol) = yC;
srData.data(:,zcol) = zC;

srDataOut = srData;
%%-------------------------------------------------------
%function srDataOut = apply3DCorrect(srData,zScaleFactor)
%%rescale Z only - for now we are not fixing xy wobble
%zcol = findSRField(srData.parInfo,'semantic','position in sample space in Z');
%
%z = srData.data(:,zcol);
%
%%apply the scale factor
%zC = z*zScaleFactor;
%
%srData.data(:,zcol) = zC;
%
%srDataOut = srData;
%-------------------------------------------------------
function [srData xyDataPH_Pix]= transformToPH(srData,srPixSize,phPixSize,tformFname,imParamPh,imParamFl)
%transform the SR data to PH nm units

m = mTrack();

xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');

%1. transform SR data to SR pixel units 
% Rapidstorm has 1 pixel offset and crops half a pixel round edges to avoid coord system confusion
xyData = srData.data(:,[xcol, ycol]);
xyData(:,1) = xyData(:,1)./srPixSize(1) + 1;
xyData(:,2) = xyData(:,2)./srPixSize(2) + 1;

%2. Account for SR coord system offset
xyData(:,1) = xyData(:,1)+imParamFl.ROI(1)-1;
xyData(:,2) = xyData(:,2)+imParamFl.ROI(2)-1;
%3. apply transfor to PH system
datatmp = load(tformFname,'dblTform_phToFl');
dblTform_phToFl= datatmp.dblTform_phToFl;
xyDataPH = m.dblTformInv(dblTform_phToFl,xyData);
%change points at image edges which may not have been sucessfully transformed to nan
xyDataPH = nanOutsidePh(xyDataPH,imParamPh.ROI);


%4. Account for PH coord system offset
xyDataPH(:,1) = xyDataPH(:,1) - (imParamPh.ROI(1)-1);
xyDataPH(:,2) = xyDataPH(:,2) - (imParamPh.ROI(2)-1);

%For debugging and visualisation
xyDataPH_Pix=xyDataPH;

%5. transform back to nanometres
xyDataPH(:,1) = xyDataPH(:,1).*phPixSize(1);
xyDataPH(:,2) = xyDataPH(:,2).*phPixSize(2);
%put the SR data structure in the new coord system
srData.data(:,xcol) = xyDataPH(:,1);
srData.data(:,ycol) = xyDataPH(:,2);

%----------------------------------------------------------------
function [srDataOut , xDriftNm,yDriftNm]= driftCorrectPH(srData,phPathPre,phPathPost,nFrameFl,phPixSize);
%linear drift correction based on drift in phase contrast images


% calculate the drift correction
phPre = imread(phPathPre);
phPost = imread(phPathPost);
[xDriftPix yDriftPix] = getDriftPH(phPre,phPost); 
xDriftNm =xDriftPix*phPixSize(1);
yDriftNm =yDriftPix*phPixSize(2);

%apply the drift correction
frameDriftX = @(fr) fr* xDriftNm/(nFrameFl-1);
frameDriftY = @(fr) fr* yDriftNm/(nFrameFl-1);

xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
fcol = findSRField(srData.parInfo,'semantic','frame number');
x = srData.data(:,xcol);
y = srData.data(:,ycol);
f = srData.data(:,fcol);

xC = x - frameDriftX(f);
yC = y - frameDriftY(f);

srDataOut = srData;
srDataOut.data(:,xcol)=xC;
srDataOut.data(:,ycol)=yC;

