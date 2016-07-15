function zPos = getZPosSimple(sigmaXY_nm, calData)

sigmaXY_um = sigmaXY_nm./1e3;
zF = calData(:,1);
sxF = calData(:,2);
syF = calData(:,3);

nPoint = numel(zF)-1;
%calculate the spline
ppZ =splinefit(zF',[sxF,syF]',nPoint,'r');
zLim = [min(zF),max(zF)];

nPts = size(sigmaXY_um,1);
zPos= zeros(nPts,1);
for ii = 1:nPts
   sx = sigmaXY_um(ii,1);
   sy = sigmaXY_um(ii,2);
   zPos(ii) = findNearestZ(sx,sy,ppZ,zLim);
end

%
%figure;plot3(sxF,syF,zF, 'ro')
%zO = zLim(1):10:zLim(2);
%v = ppval(ppZ,zO);
%sxO  = v(1,:);
%syO  = v(2,:);
%hold all;
%plot3(sxO,syO,zO,'r-');
%xlabel('sx');
%ylabel('sy');
%zlabel('z');
%------------------------------------
function z = findNearestZ(sxF,syF,ppZ,zLim)
%do a least square fit to find z

initGuessZ = mean(zLim);
lb = zLim(1);
ub = zLim(2);

f=@(z,ppZ) ppval(ppZ,z);
ydata = [sxF;syF];

options = optimset('Display','off');
z = lsqcurvefit(f,initGuessZ,ppZ,ydata,lb,ub,options);


