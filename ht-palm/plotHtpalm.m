function  plotHtpalm(htPalmSummary)

%%load all the data into a big cell
%[srInMesh htPalmConfigData]= loadAllData(htPalmSummary);
%
%srInMesh = deleteBadMesh(srInMesh);
%save('tempResult.mat');
load('tempResult.mat');



gammaVal = .7;
satVal = 0.001;
%plot kymograph
figName = [htPalmConfigData.acqFolder,filesep,htPalmConfigData.acqFileBaseName,'_kymograph.fig'];
%renderKymograph(srInMesh,figName,'Gamma',gammaVal,'Saturate',satVal);

%plot 2D movie
gammaVal = 1;
satVal = 0.001;
movieName = [htPalmConfigData.acqFolder,filesep,htPalmConfigData.acqFileBaseName,'_2dSrMovieInt_pix10sigma10'];
figFolder= [htPalmConfigData.acqFolder,filesep,htPalmConfigData.acqFileBaseName,'_2dSrMovieInt_pix10sigma10_outline'];
figName= [htPalmConfigData.acqFileBaseName,'_2dSrMovieInt_pix10sigma10_outline'];
plotHtpalm2DInt(srInMesh, movieName,gammaVal,satVal,figFolder,figName);

%plot d-z cross section
gammaVal = 1;
satVal = 0.001;
lRange = [-600 300];
movieName = [htPalmConfigData.acqFolder,filesep,htPalmConfigData.acqFileBaseName,'_2dSrMovieCrossRaw_pix10sigma10.tif'];
%plotHtpalmCrossSection(srInMesh, movieName,lRange,gammaVal,satVal);

ringStatName= [htPalmConfigData.acqFolder,filesep,htPalmConfigData.acqFileBaseName,'_ringStats.mat'];
%calculateRingRadius(srInMesh,lRange,ringStatName);

