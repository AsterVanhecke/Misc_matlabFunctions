function htPalmSummary = analyseMTS(dataFolder,dataConfig,pixSizeRescale, mTrackParam,zScaleFactor,xyWobbleFile,doFlipZ,tformName,varargin);

isTestRun = false;
is2D = false;
ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'TestRun')
    isTestRun = true;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'2DFit')
   is2D = true;
    ii=ii+1;
 else
    ii=ii+1;
 end
end


%% load the config data
[configData, fnameList] = loadConfigData(dataFolder,dataConfig,isTestRun);
fnameList = getAcqTime(fnameList,mTrackParam.tExpt_min);

configData.acqFolder = dataFolder;

configData.xyWobbleFile = xyWobbleFile;
configData.zScaleFactor = zScaleFactor;
configData.pixSize= pixSizeRescale

% figure out the names of :
%  all phPre/phPost/Sr triplets
%   **  this is fnameList
%  all SR files (single list for batch analysis)
%  all ph files (list of pairs of pre/post)

%% rescale and wobble correct the data
%TODO
display('Load data & correct abberation');
fnameList = wobbleCorrectMTS(fnameList,isTestRun,xyWobbleFile,zScaleFactor,pixSizeRescale,doFlipZ,is2D);
%drift correction
display('Correct drift');
fnameList = driftCorrectMTS(fnameList,isTestRun,mTrackParam,configData,tformName);
%% run the microbetracker analysis
%TODO
display('Manual segmentation');
fnameList = mtsSegmentation(fnameList,mTrackParam,isTestRun);

%at this point expect cell containing:
%  htPalmFiles{ii}.
%     folderName
%     phPre
%     phPost
%     srRawName
%     srLoc
%     phPreMeshName
%     srInMeshName

display('Put locs in mesh');
% combine SR localizations and microbe tracker, save the data
fnameList = putMtsLocInMeshAll(fnameList,mTrackParam,dataFolder, isTestRun);

%fnameList = phMTS_locCombine(fnameList, configData,isTestRun);
if isTestRun
  htPalmSummary= [configData.acqFolder,filesep,configData.acqFileBaseName,'_htpalmSummary.TESTRUN.mat'];
else
  htPalmSummary= [configData.acqFolder,filesep,configData.acqFileBaseName,'_htpalmSummary.mat'];
end
save(htPalmSummary,'fnameList','configData');

%-------------------------------------------------------
function [configData, fnameList] = loadConfigData(acqFolder,acqConfigFile,isTestRun);


% import the xml file
configData.config = xml_read([acqFolder,filesep(),acqConfigFile]);% using xml_read from xml_io_tools_2010_11_0
configData.acqFileBaseName = configData.config.fileBaseName_;
configData.acqMetadataFname = [acqFolder,filesep(),configData.acqFileBaseName,'_htpalmMetadata.xml'];
configData.acqMetadata = xml_read(configData.acqMetadataFname);

configData.folderNoStub = '_1';
configData.imStub = '_MMImages.ome.tif';
configData.imMetadataStub = '_MMImages_metadata.txt';

% get the list of fl images
% and the list of phPre images
if isTestRun
  nAcq=3;
else
  nAcq = configData.acqMetadata.nAcq_;
end
  nAcq
ii =1;
moreFilesExist = 1;
while ii <=nAcq & moreFilesExist
   fnameList{ii}.flName   = [configData.acqMetadata.acqMetadataList_.acqMetadata(ii).acqNameFl(1).CONTENT];
   fnameList{ii}.flFolder   = fullfile(acqFolder,[fnameList{ii}.flName,configData.folderNoStub]);

   if exist(fnameList{ii}.flFolder,'dir')
     fnameList{ii}.flImName = fullfile(fnameList{ii}.flFolder, [fnameList{ii}.flName,configData.imStub]);
     fnameList{ii}.flImMeta =  fullfile(fnameList{ii}.flFolder, [fnameList{ii}.flName,configData.imMetadataStub]);

     fnameList{ii}.phPreName =[configData.acqMetadata.acqMetadataList_.acqMetadata(ii).acqNamePhPre_];
     fnameList{ii}.phPreFolder = fullfile(acqFolder,[fnameList{ii}.phPreName,configData.folderNoStub]);
     fnameList{ii}.phPreImName = fullfile(fnameList{ii}.phPreFolder,[fnameList{ii}.phPreName,configData.imStub]);
     fnameList{ii}.phPreImMeta = fullfile(fnameList{ii}.phPreFolder,[fnameList{ii}.phPreName,configData.imMetadataStub]);

     fnameList{ii}.phPostName = [configData.acqMetadata.acqMetadataList_.acqMetadata(ii).acqNamePhPost_];
     fnameList{ii}.phPostFolder = fullfile(acqFolder,[fnameList{ii}.phPostName, configData.folderNoStub]);
     fnameList{ii}.phPostImName = fullfile(fnameList{ii}.phPostFolder,[fnameList{ii}.phPostName,configData.imStub]);
     fnameList{ii}.phPostImMeta = fullfile(fnameList{ii}.phPostFolder,[fnameList{ii}.phPostName,configData.imMetadataStub]);
     ii = ii+1;
   else %the acquisition has been terminated early
     fnameList(ii) = [];%delete that element
     moreFilesExist =0;%break out of loop, this is assumed to be the last file
   end
end

% get the pixel size for each camera
configData.imSummaryPh = readImSummaryMM(fnameList{1}.phPreImName);
configData.imSummaryFl = readImSummaryMM(fnameList{1}.flImName);
configData.pixUmPh = configData.imSummaryPh.PixelSize_um;
configData.pixUmFl = configData.imSummaryFl.PixelSize_um;
configData.roiPixPh = configData.imSummaryPh.ROI;
configData.roiPixFl = configData.imSummaryFl.ROI;

%-----------------------------------------------------------------
function fnameListOut = getAcqTime(fnameList,t0)
nAcq = numel(fnameList);
for ii = 1:nAcq
   data = loadjson(fnameList{ii}.phPreImMeta);
   fnameList{ii}.timeString = data.FrameKey_0x2D_0_0x2D_0_0x2D_0.Time;
   if ii == 1
      fnameList{ii}.tExpt_min = t0;
   else
      fnameList{ii}.tExpt_min = t0+ getExptTime(fnameList{1}.timeString,fnameList{ii}.timeString);
   end
   tExpt(ii) = fnameList{ii}.tExpt_min;
end

%reorder fnameList so analysis is always performed in time order
fnameListOut= cell(size(fnameList));
[tSort idx] =sort(tExpt) 
for ii = 1:nAcq
  fnameListOut{ii} = fnameList{idx(ii)};
end

fnameListOut = fnameList;
%-------------------------------------------------------------------
function tExpt_min = getExptTime(tStringStart,tStringCur);
n0=datenum(tStringStart,'yyyy-mm-dd HH:MM:SS');
nC=datenum(tStringCur,'yyyy-mm-dd HH:MM:SS');

tExpt_min = (nC-n0)*1440;%days to minutes

%-------------------------------------------------------------------
function fnameListOut= wobbleCorrectMTS(fnameList,isTestRun,xyWobbleFile,zScaleFactor,pixSizeRescale,doFlipZ,is2D)

%function fnameListOut = srAnalysis(fnameList,srPixSize, SNR, calFileName,isTestRun,doSkipRapidStorm,doDeletePalmIm,nDeleteSkip)

%make the fname list into a single cell list of file names
% and assign the SR output file name
if isTestRun
  nAcq=3;
else
  nAcq = numel(fnameList);
end

for ii =1:nAcq
  ii
   fnameList{ii}.srRawName = [fnameList{ii}.flImName(1:end-3),'txt']; %NB - if rapidstorm ever changes its naming scheme will need to change this!!
   fnameList{ii}.srCorName = [fnameList{ii}.flImName(1:end-3),'mat']; %wobble corrected and rescaled in xyz

   correctAbberation(fnameList{ii}.srRawName,fnameList{ii}.srCorName,xyWobbleFile,zScaleFactor,pixSizeRescale,doFlipZ,is2D);

end

fnameListOut = fnameList;
%--------------------------------------------
function correctAbberation(srRawName,srCorName,xyWobbleFile,zScaleFactor,pixSizeRescale,doFlipZ,is2D);

% load the SR analysis
srData = readRapidStorm(srRawName);
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
fcol = findSRField(srData.parInfo,'semantic','frame number');
photCol =  findSRField(srData.parInfo,'semantic','emission strength');

x = srData.data(:,xcol);
y = srData.data(:,ycol);

if ~is2D
  zcol = findSRField(srData.parInfo,'semantic','position in sample space in Z');

  z = srData.data(:,zcol);
  %wobble correction and Z-rescale
  [xC,yC,zC] = correct3D(x,y,z,xyWobbleFile,zScaleFactor);

  if doFlipZ
    zC = -1*zC;
  end

  srData.data(:,zcol) = zC;
else
  xC=x;
  yC=y;
end

%rescale the pixel size
xC = xC*pixSizeRescale(1);
yC=  yC*pixSizeRescale(2);

srData.data(:,xcol) = xC; 
srData.data(:,ycol) = yC;

save(srCorName,'srData');
%--------------------------------------------
 function fnameListOut = mtsSegmentation(fnameList,param,isTestRun);

%for each file
%make the fname list into a single cell list of file names
% and assign the SR output file name
if isTestRun
  nAcq=3;
else
  nAcq = numel(fnameList);
end

for ii =1:nAcq
   fnameList{ii}.srCorName = [fnameList{ii}.flImName(1:end-3),'mat']; 
   fnameList{ii}.srSegRoiName = [fnameList{ii}.flImName(1:end-3),'_segRoi.mat']; 
   fnameList{ii}.meshName = [fnameList{ii}.flImName(1:end-3),'_mesh.mat']; 

end

h1=figure;
h2 = figure;
%semi-manually segment the image (record the operations)
 for ii =1:nAcq
  ii
   fnameList{ii}.tExpt_min
   fnameList{ii}.meshName
%   Generate the SR image for segmentation
   srIm = makeSegIm(fnameList{ii}.srCorName,param.srImPixSize,param.srImLim,param.srBlurSigma);
   srImSat = imadjust(srIm,stretchlim(srIm),[0 1]);
   figure(h1);imshow(srImSat);
   [srImSeg recordSteps] = semiAutoSeg(srIm,param,h2);
  save(fnameList{ii}.srSegRoiName,'srIm','srImSeg','recordSteps');
 end

pixSize = param.pixSize;
h3= figure;
for ii =1:nAcq
  ii
  try
    load(fnameList{ii}.srSegRoiName,'srIm','srImSeg');
    cellList = getMeshAlg1(srIm,srImSeg,param);
    save(fnameList{ii}.meshName,'cellList','pixSize') ;
    %figure(h3);
    %dispcellall(cellList,srIm);
  catch ME
   cellList = cell(1);
   cellList{1} = cell(1);
    save(fnameList{ii}.meshName,'cellList','pixSize') ;
   fprintf('\n**************************************************\n');
   fprintf(  '******* File skipped for following reason: *******\n');
   fprintf('%s\n',getReport(ME));
   fprintf(  '**************************************************\n');
  end
end

%save the cells
fnameListOut = fnameList;
 %end
%------------------------------------------------------------------------------
function srIm = makeSegIm(srCorName,srImPixSize,srImLim,srBlurSigma);

% load the SR analysis
load(srCorName,'srData');
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
x = srData.data(:,xcol);
y = srData.data(:,ycol);

srIm = stormHist2d(x,y,srImPixSize,'XLim',srImLim{1},'YLim',srImLim{2},'BlurSigma',srBlurSigma,'Normalize',0.0);
srIm = im2uint8(srIm);

%------------------------------------------------------------------------------
function fnameList = putMtsLocInMeshAll(fnameList,mTrackParam,dataFolder,isTestRun);

%for each file
%make the fname list into a single cell list of file names
% and assign the SR output file name
if isTestRun
  nAcq=3;
else
  nAcq = numel(fnameList);
end

for ii =1:nAcq
  ii
   fnameList{ii}.srInMeshName =[dataFolder, filesep(),fnameList{ii}.flName,'_SR_Mesh.mat'];
   putMtsLocInMesh( fnameList{ii}.meshName,fnameList{ii}.srInMeshName, fnameList{ii}.srCorName,mTrackParam);
end

%------------------------------------------------------------------------------
function putMtsLocInMesh( meshName,srInMeshName, srCorName,mTrackParam);

% load the SR analysis
load(srCorName,'srData','driftFlNm','flDriftTotal');

try 
  meshData = load(meshName);
  cellList = meshData.cellList;
  cellList= changeUnitSystem(cellList, mTrackParam.srImPixSize,'nanometres');

  % put SR points in PH meshes (use spotfinder routines)
  cellList0=cellList;
  cellList = associatePhSr(cellList, srData);
  cellList = addExptTime(cellList,mTrackParam.tExpt_min);
  cellList = addDriftRecord(cellList,driftFlNm,flDriftTotal);

  %%<DEBUG
  %%plot the image and localisations
  %xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
  %ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
  %x=srData.data(:,xcol);y=srData.data(:,ycol);
  %figure;hold all
  %plot(x,y,'r.');
  %nCell = numel(cellList{1});
  %ii=1
  %for jj =1:nCell
  %  if ~isempty(cellList{ii}{jj})
  %    mesh = cellList{ii}{jj}.mesh;
  %    plot(mesh(:,1),mesh(:,2),'g-');%left side of cell
  %    plot(mesh(:,3),mesh(:,4),'g-');%right side of cell
  %  end
  %end
  %axis equal
  %pause
  %%DEBUG>
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

%---------------------------------------------------------------------
function fnameList = driftCorrectMTS(fnameList,isTestRun,mTrackParam,configData,tformName)

phPixSizeNm =configData.pixUmPh*1e3; %CONVERTED to nm
nFrameFl = configData.config.camEmccdNumFrames_;
flPixSizeNm = mTrackParam.pixSize;

if isTestRun
  nAcq=3;
else
  nAcq = numel(fnameList);
end

for ii =1:nAcq
  ii
  srFile = fnameList{ii}.srCorName;
  phPathPre = fnameList{ii}.phPreImName;
  phPathPost = fnameList{ii}.phPostImName;
  driftCorrectFl(srFile,phPathPre,phPathPost,nFrameFl,phPixSizeNm,flPixSizeNm,tformName);
end

fnameListOut = fnameList;
%-------------------------------------------------------------
function driftCorrectFl(srCorFile,phPathPre,phPathPost,nFrameFl,phPixSizeNm,flPixSizeNm,tformName);

%linear drift correction based on drift in phase contrast images

% calculate the drift correction
phPre = imread(phPathPre);
phPost = imread(phPathPost);
[xDriftPix yDriftPix] = getDriftPH(phPre,phPost); 
%xDriftNm =xDriftPix*phPixSize(1);
%yDriftNm =yDriftPix*phPixSize(2);
%transform drift to fluorescence coordinate system

%transform the loc points back to the Fl coord system (nm)
m = mTrack();
%3. apply transfor to PH system
load(tformName);
TFORM_TYPE = 'lwm';
dblTform_flToPh= m.getDblTform(posPixFl_ROI,posPixPh_ROI,TFORM_TYPE);
%place the drift coords in centre of transform system
pos0 = mean(posPixPh,1);
pos1 = pos0+ [xDriftPix,yDriftPix];
posPH = [pos0;pos1];
posFl =  dblTformInv(dblTform_flToPh,posPH); 

driftFlPix = posFl(2,:)-posFl(1,:);
driftFlNm = driftFlPix.*flPixSizeNm;

%phDriftTotal = sqrt(xDriftPix^2+yDriftPix^2)*phPixSizeNm
flDriftTotal = sqrt(driftFlNm(1)^2 + driftFlNm(2)^2)

%figure;
%hold off;
%quiver(0,0,xDriftPix*phPixSizeNm,yDriftPix*phPixSizeNm);
%hold all;
%quiver(0,0,driftFlNm(1),driftFlNm(2));
%legend('PH','Fl');
%
%figure;
%plot(posPixPh(1:3,1)-posPixPh(1,1),posPixPh(1:3,2)-posPixPh(1,2),'x-');
%hold all;
%plot(posPixFl(1:3,1)-posPixFl(1,1),posPixFl(1:3,2)-posPixFl(1,2),'x-');
%legend('PH','Fl')

%apply the drift correction
load(srCorFile,'srData');
frameDriftX = @(fr) fr* driftFlNm(1)/(nFrameFl-1);
frameDriftY = @(fr) fr* driftFlNm(2)/(nFrameFl-1);

xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
fcol = findSRField(srData.parInfo,'semantic','frame number');
x = srData.data(:,xcol);
y = srData.data(:,ycol);
f = srData.data(:,fcol);

xC = x - frameDriftX(f);
yC = y - frameDriftY(f);

srData.data(:,xcol)=xC;
srData.data(:,ycol)=yC;

save(srCorFile,'srData','driftFlNm','flDriftTotal');

%--------------------------------------------------------
function cellListOut = addDriftRecord(cellList,driftFlNm,flDriftTotal)
cellListOut = cellList;

for ii = 1:numel(cellListOut{1})
   if ~isempty(cellListOut{1}{ii})
      cellListOut{1}{ii}.driftTotal = flDriftTotal;
      cellListOut{1}{ii}.drift= driftFlNm;
      cellListOut{1}{ii}.isDriftCor = true;
   end
end


