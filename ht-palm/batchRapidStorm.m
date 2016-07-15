function batchRapidStorm(fname, varargin)
%Example of working arguments
% rsarg={'--FileType','TIFF','','--PixelSizeInNM','117,126','--FitJudgingMethod','SquareRootRatio','--SNR','50','--ThreeD','Spline3D','--ZCalibration',calFileName,'','--ChooseT','Table','--ChooseT','Image','--ColourScheme','ByCoordinate','--HueCoordinate','PositionZ','--AutoTerminate','--Run'};
% RS command line is mental.
% Arguments change depending on previous arguments, cannot get a full list!
% Need to run --help each time you add a new argument to see if any new ones have appeared
% --AutoTerminate only works if you add '--ChooseT','Image' as arguments - otherwise execution just hangs
% NB: if file paths (fname, or calFileName) have spaces in them, they must be wrapped in double quotes, 
%  eg '"my awesome image.tif"'

narg =  numel(varargin);
rsargStr = '';
for ii = 1:narg
   rsargStr = [rsargStr,' ',varargin{ii}];
end

if ~iscell(fname)
   rsargStrFull = [' --InputFile ',fname,rsargStr];
   [status, result] =system(['rapidstorm.exe',rsargStrFull])
else
   for ii = 1:numel(fname)
      fCur = fname{ii};
      rsargStrFull = [' --InputFile ',fCur,rsargStr];
      [status, result] =system(['rapidstorm.exe',rsargStrFull]);
   end
end

