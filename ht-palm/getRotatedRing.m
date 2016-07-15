function [x_th y_rad meshRot] = getRotatedRing(ringRange,x,y,meshAligned,pixSize,steplength,srInMesh)
%TODO Export the rotated mesh as well

%find the mesh points within the ring range
centreLine = [mean(meshAligned(:,[1,3]),2),mean(meshAligned(:,[2,4]),2)];
ringRangeIdx = getClosest(ringRange,centreLine(:,1));
nStep = size(centreLine,1);
%ends of the cells cause problems
%so dont just take the last 2 points
%and skip the end point
if ringRangeIdx(1) ==1
  ringRangeIdx(1) = 2;
  ringRangeIdx(2) = 4;
elseif ringRangeIdx(2) == nStep
  ringRangeIdx(1) = nStep-3;
  ringRangeIdx(2) = nStep-1;
elseif ringRangeIdx(2) == ringRangeIdx(1)
  ringRangeIdx(2) = ringRangeIdx(2)+1;
end

%Get the angle of the mesh centreline
%Get the raw (unaligned mesh)
meshRaw = srInMesh.mesh;
cLineRaw =  [mean(meshRaw(:,[1,3]),2),mean(meshRaw(:,[2,4]),2)];
ring_cLine_start = cLineRaw(ringRangeIdx(1),:);
ring_cLine_end   = cLineRaw(ringRangeIdx(2),:);

cLineCtr=mean([ring_cLine_start;ring_cLine_end],1)';
dX = [ring_cLine_end-ring_cLine_start];

theta = atan2(dX(2),dX(1));

%rotate x & y by that angle
R = [cos(-theta), -sin(-theta);sin(-theta),cos(-theta)];
X = R*[x';y'];
cLineCtrRot = R*cLineCtr;

%supply the rotated mesh as well
meshTop = meshRaw(:,1:2)';
meshBottom= meshRaw(:,3:4)';
meshTop = R*meshTop;
meshBottom = R*meshBottom;
meshRot = [meshTop', meshBottom'];


%remove the offset
x_th = X(1,:)'-cLineCtrRot(1);
y_rad = X(2,:)'-cLineCtrRot(2);
meshRot(:,1) = meshRot(:,1)-cLineCtrRot(1);
meshRot(:,3) = meshRot(:,3)-cLineCtrRot(1);
meshRot(:,2) = meshRot(:,2)-cLineCtrRot(2);
meshRot(:,4) = meshRot(:,4)-cLineCtrRot(2);

%%DEBUG
%hold off;
%plot(meshRaw(:,1),meshRaw(:,2),'k-');
%hold on;
%plot(meshRaw(:,3),meshRaw(:,4),'k-');
%plot(cLineRaw(:,1),cLineRaw(:,2),'k-');
%ringLine = [ring_cLine_start;ring_cLine_end];
%plot(ringLine(:,1),ringLine(:,2),'r-');
%axis equal
%theta/pi*180
%pause
%---------------------------------------
function idx = getClosest(valToFind,data)

for ii = 1:numel(valToFind)
 [c idx(ii)] = min(abs(data-valToFind(ii)));
end
