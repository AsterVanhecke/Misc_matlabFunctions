function [xC,yC,zC] = correctWobble(x,y,z,varargin)
% CORRECTWOBBLE Correct z-dependent axial (x,y) localization error
%
% SYNTAX 
%  [xC,yC] = correctWobble(x,y,z,'WobbleMatrix',wobbleData) 
%  [xC,yC] = correctWobble(x,y,z,'WobbleMatrix',wobbleData)
%  [xC,yC] = correctWobble(x,y,z,'WobbleFile',wobbleFileName) 
%  [xC,yC,zC] = correctWobble(... ,'ZScaleFactor',zScaleFactor) 
%
% DESCRIPTON
%  [xC,yC] = correctWobble(x,y,z,'WobbleMatrix',wobbleData) returns corrected
%     Returns corrected values (xC,yC) for localizations (x,y,z), given a correction matrix
%     'wobbleData', containing a list of axial shifts as a function of Z, ie 
%     wobbleData = [Z1, xShift1, yShift1;...
%                   Z2, xShift2, yShift2;... etc]
%
%     The wobbleData (or wobbleFile) should be generated using the wobbleCalibration()
%     function
% [xC,yC] = correctWobble(x,y,z,'WobbleFile',wobbleFileName) loads a file named
%     wobbleFileName containing the matrix wobbleData, as defined above.
%     (Optionally, the xShift/ yShift columns can be columns 4,5 of the matrix,
%     instead of columns 2,3 of the matrix for legacy compatibility)
% [xC,yC,zC] = correctWobble(x,y,z,... ,'ZScaleFactor',zScaleFactor) optionally 
%     performs linear z-axis correction of spherical abberation, by a factor of
%     zScaleFactor (Huang et al, Nat. Methods 2008 suggest zScaleFactor=0.72)
%
% EXAMPLE
%  Test wobble correction on the test data supplied with these functions
%  The dataset below ('bead 0.1um -1.5to1.5 20nmstep small FOV test correction.txt')
%  shows a fluorescent bead (Invitrogen 0.1um TetraSpeks),
%  stepped in Z with a piezo stage and 3D localized with RapidSTORM (in X,Y and Z). 
%
%  Define the file names (this is the test data supplied with the 
%  wobble correction functions):
%     fname = 'test data\bead 0.1um -1.5to1.5 20nmstep small FOV test correction.txt';
%     wobbleName = 'test data\150309 cal bead 0.1um -1.5to1.5 20nmstep.range-725To892.xyWobble.txt';
%  Load the data
%      a=importdata(fname);a=a.data;
%      x=a(:,2);y=a(:,3);z=a(:,4);
%  Run the correction
%      [xC,yC] = correctWobble(x,y,z,'WobbleFile',wobbleName);
%  Compare the corrected and uncorrected data via scatter plot 
%      figure;hold all;
%      plot(x,y,'r.');
%      plot(xC,yC,'kx');
%      axis equal
%      legend('Raw','Corrected');
%      xlabel('X (nm)');
%      ylabel('Y (nm)');
%       
% This software is released under the GPL v3 (see license file 'gpl.txt'). It is provided AS-IS and no
% warranty is given.
%
% Author: Seamus Holden
% Last update: April 2015

narg = numel(varargin);
zScaleFactor = 1;%DEFAULT - no correction
wobbleData=[];
ii=1;
while ii<=narg
   if strcmp(varargin{ii},'WobbleMatrix')
      wobbleData = varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'WobbleFile')
      wobbleFile = varargin{ii+1};
      wobbleData = importdata(wobbleFile);
      ii = ii+2;
   elseif strcmp(varargin{ii},'ZScaleFactor')
      zScaleFactor = varargin{ii+1};
      ii = ii+2;
   else 
      ii=ii+1;
   end
end

if isempty(wobbleData)
   error('No wobble correction data supplied');
end

if size(wobbleData,2) == 3
   zF = wobbleData(:,1);
   xF = wobbleData(:,2);
   yF = wobbleData(:,3);
elseif size(wobbleData,2) == 5
   zF = wobbleData(:,1);
   xF = wobbleData(:,4);
   yF = wobbleData(:,5);
else
   error('Wobble data has wrong number of columns');
end

%apply measured correction for xy wobble, and rescale z due to spherical abberation

nPoint = numel(zF)-1;
%calculate the spline
ppX =splinefit(zF,xF,nPoint,'r');
ppY =splinefit(zF,yF,nPoint,'r');

%apply the wobble correction to values within the z limits
zlim = [min(zF), max(zF)];
inZ = z>zlim(1) & z<zlim(2);
zLow = z<zlim(1);
zHigh= z>zlim(2);

%apply the wobble correction
xC(inZ) = x(inZ) - ppval(ppX,z(inZ));
yC(inZ) = y(inZ) - ppval(ppY,z(inZ));
xC(zLow) = x(zLow) - ppval(ppX,zlim(1));
yC(zLow) = y(zLow) - ppval(ppY,zlim(1));
xC(zHigh) = x(zHigh) - ppval(ppX,zlim(2));
yC(zHigh) = y(zHigh) - ppval(ppY,zlim(2));

%apply the scale factor
zC = z*zScaleFactor;
