function [x,y,mesh] = alignCM(x,y,mesh)
%recentre and flip the data along CofM

%recentre the data around the middle of the cell
xC = mean([mesh(:,1);mesh(:,3)]);

mesh(:,1) = mesh(:,1)-xC;
mesh(:,3) = mesh(:,3)-xC;
x = x-xC;

%flip so that CM is [top-left] - NB for publication, need to make sure handedness is maintained in any 3D data!
xCM = mean(x);
if xCM>0
  mesh(:,1) = -mesh(:,1);
  mesh(:,3) = -mesh(:,3);
  x = -x;
  assert(mean(x)<0,'Error, axis not flipped sucessfully');
end

yCM = mean(y);
if yCM>0
  mesh(:,2) = -mesh(:,2);
  mesh(:,4) = -mesh(:,4);
  y = -y;
  assert(mean(y)<0,'Error, axis not flipped sucessfully');
end 


