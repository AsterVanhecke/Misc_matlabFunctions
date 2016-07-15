%-----------------------------------------
%ALL CELL DETECTION/ MESH GENERATION PARAMETERS
%image parameters
param.invertimage=1;
param.erodeNum = 0;
%segmentation parameters
param.edgemode='sauvola';
param.localNeighbourhood=[15 15];
param.localThresh = 0.2;
param.smoothNum= 10;
param.areaMin = 1000; 


% Shape parameters
param.interpoutline =0;
param.interpWeights =[0.5 0.5];
param.interpSigma =[5 5];
param.thresFactorF = 1;

% Mesh parameters;
param.getmesh = 1;
param.meshStep = 10;
param.meshTolerance = 0.01;
param.meshWidth=50;
param.fsmooth=Inf;
param.fmeshstep =10;

% Limits on poition / number;
param.roiBorder =5;
param.noCellBorder = 5;
param.imRoi = [1 1 Inf Inf];
%-----------------------------------------

%Super-res paramaters
param.srImPixSize = 10;
param.srImLim{1} = [0 (128-1)*120]; 
param.srImLim{2} = [0 (126)*120]; 
pixSizeRescale =[1 1];%CORRECT FOR WRONG RS PIXEL SIZE
param.pixSize = [120 120];
param.srBlurSigma=30;%MAYBE CHANGE ME
param.tExpt_min = 56;%ALWAYS CHANGE ME !!!!!!!!!!!!!!!!!!
doFlipZ = true;
xyWobbleFile ='G:\cameraSettings\astigmatic 3d calibration\150309 cal bead 0.1um -1.5to1.5 20nmstep.zCal.range-725To892.xyWobble.txt';
zScaleFactor=0.72;
dataFolder='G:\Anna\mts2d_4';%CHANGE ME
dataConfig ='mts-1pcA-smallCh-acqRm-dend2-3v-htpalm-10ms-al11-100-200_t0-56min_config.xml';%CHANGE ME

%Generate the 2ch tform - for drift correction
%first test the alignment data
beadFolder = 'G:\seamus-highThroughput\150622_mts_2d_repeats\2chCal_1';%CHANGE ME
beadConfigFile = '2chCal_config.xml';%CHANGE ME
saveFolder = beadFolder;
tformSaveName = [beadFolder,filesep(),'beadCal561nmPH_TFORM.mat'];
showFits=false;

%2 CHANNEL CALIBRATION
% autoAlignChannel(beadFolder,beadConfigFile,saveFolder, tformSaveName,showFits);

%MTS ANALYSIS
% htPalmSummary = analyseMTS(dataFolder,dataConfig,pixSizeRescale, param,zScaleFactor,xyWobbleFile,doFlipZ,tformSaveName,'2DFit');

 htPalmSummary='G:\Anna\mts2d_4\mts-1pcA-smallCh-acqRm-dend2-3v-htpalm-10ms-al11-100-200_t0-56min_htpalmSummary.mat';

%load all the data into a big cell
[srInMeshAll htPalmConfigDataAll]= loadAllData(htPalmSummary);
% save('tmpResult150615Mts150611_1.mat');
% load('tmpResult150615Mts150611_1.mat');

%filter the cells based on manual classification
% blurSigma=10;%MIGHT WANT TO ALTER THIS
% manuallyClassifyDivision(srInMeshAll, htPalmSummary,param.tExpt_min,blurSigma,param.srImLim);
 manClassifyFile=[htPalmSummary(1:end-4),'.manualClassify.mat'];

%% Extract radii for all cells
lStep =12.5;
lWidth =100;
histStep = 50;
showFit=true;
% Diameters of the bacteria as a function of their lengths and all the
% other info regarding each bacteria (the variable is a cell)
[srInMeshAll_fineStep]=getCellWidthAll(srInMeshAll,lStep,lWidth,histStep,showFit);
% showFit=false;
% [srInMeshAll_fineStep] = getMinMaxDiamAll(srInMeshAll_fineStep,showFit);
% %pltMTS_diameter(srInMeshAll_fineStep)
% save('tmpResult150615Mts150611_lWidth100.mat','srInMeshAll_fineStep');
%
% lStep =12.5;
% lWidth =50;
% histStep = 50;
% showFit=true;
% [srInMeshAll_fineStep]=getCellWidthAll(srInMeshAll,lStep,lWidth,histStep,showFit);
% showFit=false;
% [srInMeshAll_fineStep] = getMinMaxDiamAll(srInMeshAll_fineStep,showFit);
% pltMTS_diameter(srInMeshAll_fineStep)
% save('tmpResult150615Mts150611_lWidth50.mat','srInMeshAll_fineStep');
%
% lStep =12.5;
% lWidth =25;
% histStep = 50;
% showFit=true;
% [srInMeshAll_fineStep]=getCellWidthAll(srInMeshAll,lStep,lWidth,histStep,showFit);
% showFit=false;
% [srInMeshAll_fineStep] = getMinMaxDiamAll(srInMeshAll_fineStep,showFit);
% %pltMTS_diameter(srInMeshAll_fineStep)
% save('tmpResult150615Mts150611_lWidth25.mat','srInMeshAll_fineStep');
% %
load('tmpResult150615Mts150611_lWidth25.mat','srInMeshAll_fineStep');
load(manClassifyFile,'manClass');

tMax = 600;
tMax_filter=200;
nSigma = 1.5;

[srInMeshFilt cellOk]= deleteBadLength(srInMeshAll_fineStep,tMax_filter,nSigma,'OutputSrData','FilterLengthTime','ManualClassification',manClass,'MaxTime',tMax);
[t,w,d,cellData]=pltMTS_diameter(srInMeshAll_fineStep,'ManualClassification',manClass,'TMax',tMax,'FilterList',cellOk);
% save('tmpResult150615Mts150611_lStep12.5lWidth25_nSig1.5.mat','t','w','d');

