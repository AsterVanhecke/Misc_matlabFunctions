%first test the alignment data
beadFolder = 'G:\seamus\development\htpalmDev\130123_ftsZHtPalm\beadCal';
beadConfigFile = 'beadCal561nmPH_config.xml';
saveFolder = beadFolder;
tformSaveName = [beadFolder,filesep(),'beadCal561nmPH_TFORM.mat'];
showFits=false;

%autoAlignChannel(beadFolder,beadConfigFile,saveFolder, tformSaveName,showFits);

%now analyse the ht data

dataFolder = 'G:\seamus\development\htpalmDev\130123_ftsZHtPalm\ftsZ_htpalm_10ms_2' ;
dataConfig = 'ftsZhtpalm_10ms_config.xml';
srPixSize = [117,126];
srSNR = 50;
calFileName= '"G:\seamus\development\htpalmDev\astigmatic 3d calibration\3Df1000_zCal_evolve128_121204.txt"';
mTrackParamFile = 'G:\seamus\development\htpalmDev\htpalm-analysis\alg4_130204.set';
htPalmSummary = analyseHtpalm(dataFolder,dataConfig,tformSaveName,srPixSize, srSNR, calFileName,mTrackParamFile);
%save('tempresult.mat');
%load('tempresult.mat');

%plotHtpalm(htPalmSummary);
