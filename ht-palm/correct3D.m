function [xC,yC,zC] = correct3D(x,y,z,wobbleFile,zScaleFactor)
%apply measured correction for xy wobble, and rescale z due to spherical abberation

if ~isempty(wobbleFile)% isempty indicates to skip wobble correction
   wobbleData = importdata(wobbleFile);

   zF = wobbleData(:,1);
   xF = wobbleData(:,4);
   yF = wobbleData(:,5);

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
else
   xC = x;
   yC = y;
end

%apply the scale factor
zC = z*zScaleFactor;

%hold off;
%plot(x,y,'k.');
%hold all;
%plot(xC,yC,'r.');
%axis equal;
%pause
