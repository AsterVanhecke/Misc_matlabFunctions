%imName = {'"G:\seamus\development\htpalmDev\dataTesting\ftsZhtpalm_10ms_FOV145_Acq0_flCh0_1\ftsZhtpalm_10ms_FOV145_Acq0_flCh0_MMImages.ome.tif"',};
imName = {'"G:\seamus\development\htpalmDev\dataTesting\ftsZhtpalm_10ms_FOV145_Acq0_flCh0_1\ftsZhtpalm_10ms_FOV145_Acq0_flCh0_MMImages.ome.tif"','"G:\seamus\development\htpalmDev\dataTesting\ftsZhtpalm_10ms_FOV145_Acq0_flCh0_1\ftsZhtpalm_10ms_FOV145_Acq0_flCh0_MMImages.ome.tif"'};
calFileName= '"G:\seamus\development\htpalmDev\astigmatic 3d calibration\3Df1000_zCal_evolve128_121204.txt"';

rsarg={'--FileType','TIFF','','--PixelSizeInNM','117,126','--FitJudgingMethod','SquareRootRatio','--SNR','50','--ThreeD','Spline3D','--ZCalibration',calFileName,'','--ChooseT','Table','--ChooseT','Image','--ColourScheme','ByCoordinate','--HueCoordinate','PositionZ','--AutoTerminate','--Run'};

batchRapidStorm(imName, rsarg{:});
