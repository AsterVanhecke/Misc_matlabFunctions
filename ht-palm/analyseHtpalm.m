function htPalmSummary = analyseHtpalm(dataFolder,dataConfig,tformName,srPixSize, srSNR, srCalFile, mTrackParamFile,zScaleFactor,xyWobbleFile,finalRoiBorder,varargin);

isTestRun = false;
doRapidStormOnly=false;
doSkipRapidStorm=false;%ie if you've already localized
doDeletePalmIm=false;
nDeleteSkip = 30;%Keep each of the nth raw image files if we do delete
is2D = false;
fwhm=NaN;

ii = 1;
while ii <= numel(varargin)
 if strcmp(varargin{ii},'TestRun')
    isTestRun = true;
    ii = ii + 1;
 elseif strcmp(varargin{ii},'RapidStormOnly')
   doRapidStormOnly = true;
   ii = ii + 1;
 elseif strcmp(varargin{ii},'SkipRapidStorm')
   doSkipRapidStorm= true;
   ii = ii + 1;
 elseif strcmp(varargin{ii},'DeletePalm')
   doDeletePalmIm = true;
   ii = ii + 1;
 elseif strcmp(varargin{ii},'NDeleteSkip')
   nDeleteSkip = varargin{ii+1};
   ii = ii + 2;
 elseif strcmp(varargin{ii},'2DFit')
   is2D = true;
   fwhm= varargin{ii+1};
   ii = ii + 2;
 else
   ii=ii+1;
 end
end
%
%% load the config data
[configData, fnameList] = loadConfigData(dataFolder,dataConfig,isTestRun);
fnameList = getAcqTime(fnameList);
configData.srPixSize= [srPixSize] ;% Have to add this manually since MM only allows for square pixels, we effectively have rectangular pixels!
configData.srSNR = srSNR;
configData.srCalFile = srCalFile;
configData.mTrackParamFile = mTrackParamFile;
imSummaryPh = readImSummaryMM(fnameList{1}.phPreImName);
imSummaryFl = readImSummaryMM(fnameList{1}.flImName);
configData.roiPixPh = imSummaryPh.ROI;
configData.roiPixFl = imSummaryFl.ROI;

if ~doRapidStormOnly
  configData.tformName = tformName;
  configData.units='nanometers';
  configData.acqFolder = dataFolder;

  configData.xyWobbleFile = xyWobbleFile;
  configData.zScaleFactor = zScaleFactor;

  % figure out the names of :
  %  all phPre/phPost/Sr triplets
  %   **  this is fnameList
  %  all SR files (single list for batch analysis)
  %  all ph files (list of pairs of pre/post)

  % run rapidstorm
  fnameList = srAnalysis(fnameList,configData.srPixSize, configData.srSNR, configData.srCalFile,isTestRun,doSkipRapidStorm,doDeletePalmIm,nDeleteSkip,is2D,fwhm);

  %% run the microbetracker analysis
  fnameList = mTrackAnalysis(fnameList,configData.mTrackParamFile,finalRoiBorder,isTestRun);

  %at this point expect cell containing:
  %  htPalmFiles{ii}.
  %     folderName
  %     phPre
  %     phPost
  %     srRawName
  %     srLoc
  %     phPreMeshName
  %     srInMeshName

  % combine SR localizations and microbe tracker, save the data

  fnameList = phFlCombineAll(fnameList, configData,isTestRun);
  if isTestRun
    htPalmSummary= [configData.acqFolder,filesep,configData.acqFileBaseName,'_htpalmSummary.TESTRUN.mat'];
  else
    htPalmSummary= [configData.acqFolder,filesep,configData.acqFileBaseName,'_htpalmSummary.mat'];
  end
  save(htPalmSummary,'fnameList','configData');
  %save('testRoiMod150209.mat','fnameList','configData');
else
  % run rapidstorm
  fnameList = srAnalysis(fnameList,configData.srPixSize, configData.srSNR, configData.srCalFile,isTestRun,doSkipRapidStorm,doDeletePalmIm,nDeleteSkip,is2D,fwhm);
  htPalmSummary='';
end

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
  nAcq=1;
else
  nAcq = configData.acqMetadata.nAcq_;
end
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

%-------------------------------------------------------
function fnameListOut = srAnalysis(fnameList,srPixSize, SNR, calFileName,isTestRun,doSkipRapidStorm,doDeletePalmIm,nDeleteSkip,is2D,fwhm)

spotSearchEagerness = 1;
darkCurrent = 0;
countsPerPhoton= 1;
smoothingMaskSize=3;

%make the fname list into a single cell list of file names
% and assign the SR output file name
if isTestRun
  nAcq=1;
else
  nAcq = numel(fnameList);
end

srInputList = {};
for ii =1:nAcq
   srInputList = {srInputList{:}, ['"',fnameList{ii}.flImName,'"']};
   fnameList{ii}.srRawName = [fnameList{ii}.flImName(1:end-3),'txt']; %NB - if rapidstorm ever changes its naming scheme will need to change this!!
end

fnameListOut = fnameList;

pixSizeStr = [num2str(srPixSize(1)),',',num2str(srPixSize(2))];
fwhmStr = [num2str(fwhm),',',num2str(fwhm)];
if ~is2D
  %these are the defaults for 3D
  rsarg={'--FileType','TIFF','','--PixelSizeInNM',pixSizeStr,'--FitJudgingMethod','SquareRootRatio','--SNR',num2str(SNR),'--ThreeD','Spline3D','--ZCalibration', ['"',calFileName,'"'],'','--ChooseT','Table','--ChooseT','Image','--ColourScheme','ByCoordinate','--HueCoordinate','PositionZ','--AutoTerminate','--Motivation',num2str(spotSearchEagerness),'--DarkCurrent',num2str(darkCurrent),'--CountsPerPhoton',num2str(countsPerPhoton),'--MLEFitting','--SmoothingMaskSize',num2str(smoothingMaskSize),'--Run'};
else
  rsarg={'--FileType','TIFF','','--PixelSizeInNM',pixSizeStr,'--FitJudgingMethod','SquareRootRatio','--SNR',num2str(SNR),'--ThreeD','No3D','--PSF',fwhmStr,'--ChooseT','Table','--ChooseT','Image','--ColourScheme','ByCoordinate','--HueCoordinate','PositionZ','--AutoTerminate','--Motivation',num2str(spotSearchEagerness),'--DarkCurrent',num2str(darkCurrent),'--CountsPerPhoton',num2str(countsPerPhoton),'--MLEFitting','--SmoothingMaskSize',num2str(smoothingMaskSize),'--Run'};
end

if ~doSkipRapidStorm
  batchRapidStorm(srInputList, rsarg{:});
end

if doDeletePalmIm
  deletePalm(srInputList,nDeleteSkip);
end

%-------------------------------------------------------
function fnameListOut = mTrackAnalysis(fnameList,mTrackParamFile,finalRoiBorder)
%only analyse the phPre's - ph post just for drift correction (for now)


if isTestRun
  nAcq=1;
else
  nAcq = numel(fnameList);
end

for ii =1:nAcq
   fnameList{ii}.phPreMeshName = [fnameList{ii}.phPreFolder,filesep(),fnameList{ii}.phPreName,'_mesh.mat'];
   %doMTrack(fnameList{ii}.phPreFolder,fnameList{ii}.phPreMeshName,mTrackParamFile);

   %HACK TO MODIFY THE BOX SZ
   %WITHOUT CHANGING THE ANALYSIS
   fixBoxSize(fnameList{ii}.phPreFolder,fnameList{ii}.phPreMeshName,mTrackParamFile,finalRoiBorder);
end

fnameListOut = fnameList;
%--------------------------------------------------------
%--------------------------------------------------------
function doMTrack(phFolder,saveName,algFname);

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

%--------------------------------------------------------
%--------------------------------------------------------
function fnameListOut = phFlCombineAll(fnameList, configData)

tformFname = configData.tformName;
srPixSizeNm = configData.srPixSize;
phPixSizeNm = [ configData.pixUmPh,  configData.pixUmPh]*1e3; %CONVERTED to nm
units = configData.units;
roiPixPh = configData.roiPixPh;
roiPixFl = configData.roiPixFl;
nFrameFl = configData.config.camEmccdNumFrames_;


if isTestRun
  nAcq=1;
else
  nAcq = numel(fnameList);
end

for ii = 1:nAcq
   ii
   fnameList{ii}.srInMeshName =[configData.acqFolder, filesep(),fnameList{ii}.flName,'_SR_Mesh.mat'];
   meshFile = fnameList{ii}.phPreMeshName;
   srFile = fnameList{ii}.srRawName;
   phPathPre = fnameList{ii}.phPreImName;
   phPathPost = fnameList{ii}.phPostImName;
   xyWobbleFile = configData.xyWobbleFile;
   zScaleFactor = configData.zScaleFactor;
   tExpt_min = fnameList{ii}.tExpt_min;
   [xDriftNm, yDriftNm] = phFlCombine3(srFile,meshFile,fnameList{ii}.srInMeshName, tformFname, srPixSizeNm, phPixSizeNm, roiPixFl, roiPixPh,units,phPathPre,phPathPost,zScaleFactor,nFrameFl,tExpt_min,xyWobbleFile);
   fnameList{ii}.xDriftNm = xDriftNm;
   fnameList{ii}.yDriftNm = yDriftNm;
end
fnameListOut = fnameList;

%-----------------------------------------------------------------
function fnameListOut = getAcqTime(fnameList)
nAcq = numel(fnameList);
for ii = 1:nAcq
   data = loadjson(fnameList{ii}.phPreImMeta);
   fnameList{ii}.timeString = data.FrameKey_0x2D_0_0x2D_0_0x2D_0.Time;
   if ii == 1
      fnameList{ii}.tExpt_min = 0;
   else
      fnameList{ii}.tExpt_min = getExptTime(fnameList{1}.timeString,fnameList{ii}.timeString);
   end
end
fnameListOut = fnameList;
%-------------------------------------------------------------------
function tExpt_min = getExptTime(tStringStart,tStringCur);
n0=datenum(tStringStart,'yyyy-mm-dd HH:MM:SS');
nC=datenum(tStringCur,'yyyy-mm-dd HH:MM:SS');

tExpt_min = (nC-n0)*1440;%days to minutes

%-------------------------------------------------------------------
function fixBoxSize(phPreFolder,phPreMeshName,mTrackParamFile,finalRoiBorder)

load(phPreMeshName);
roiBorderOld =p.roiBorder;

deltaRoi = roiBorderOld-finalRoiBorder;

for ii = 1:numel(cellList{1})
  if ~isempty(cellList{1}{ii})
    mesh = cellList{1}{ii}.mesh;
    boxSmall= getBox(mesh,finalRoiBorder);
    cellList{1}{ii}.boxMeshCreation=cellList{1}{ii}.box;
    cellList{1}{ii}.box=boxSmall;

    %ii
    %hold off
    %plot(mesh(:,1),mesh(:,2))%left side of cell
    %hold all
    %plot(mesh(:,3),mesh(:,4));%right side of cell
    %plot(mesh(:,[1 3])',mesh(:,[2 4])')%segments
    %plotBox(boxSmall,'r');
    %axis equal
    %pause
  end
end
save(phPreMeshName);


%-------------------------------------------------------------------
function plotBox(box,lineColor)
if ~exist('lineColor','var');
  lineColor='k';
end
x0 = box(1);
y0 = box(2);
x1 = box(1)+box(3)-1;
y1 = box(2)+box(4)-1;
plot([x0,x1],[y0,y0] ,'-','Color',lineColor);
plot([x0,x1],[y1,y1] ,'-','Color',lineColor);
plot([x0,x0],[y0,y1] ,'-','Color',lineColor);
plot([x1,x1],[y0,y1] ,'-','Color',lineColor);
%-------------------------------------------------------------------
function  box= getBox(mesh,finalRoiBorder);
meshX=[mesh(:,1);mesh(:,3)];
meshY=[mesh(:,2);mesh(:,4)];
%top left: x0,y0; bottom right: x1,y1
x0 = min(meshX);
x1 = max(meshX);
y0 = min(meshY);
y1 = max(meshY);

box = round([x0-finalRoiBorder,y0-finalRoiBorder,x1-x0+2*finalRoiBorder,y1-y0+2*finalRoiBorder]);

%-------------------------------------------------------------------
function deletePalm(fname,nDeleteSkip)

for ii = 1:numel(fname)
  if mod(ii-1,nDeleteSkip)~=0
    f = fname{ii};
    f=f(2:end-1); %the rapidstorm arg is wrapped in quotes - remove these
    f_meta = [f(1:end-8),'_metadata.txt'];%delete the metadata too;
    delete(f);
    delete(f_meta);
  end
end

