acqFolder = 'G:\seamus\development\htpalmDev\121204_CC_Dendra2FtsZ_GrowthTest_and_BeadCal\2chBeadCal_0.5umBead_1';
acqConfigFile = 'beadCal_config.xml';
saveFolder = acqFolder;
tformSaveName = [acqFolder,filesep(),'2chBeadCal_0.5umBead_1_TFORM.mat'];
showFits=true;

autoAlignChannel(acqFolder,acqConfigFile,saveFolder, tformSaveName,showFits);
