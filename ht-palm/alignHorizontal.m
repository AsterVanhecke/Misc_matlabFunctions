function [xC,yC,meshC, theta] = alignHorizontal(x,y,mesh)
%align the data along the horizontal
%theta in radians

%align by first and last points of the mesh
X0 = [mesh(1,1), mesh(1,2)];
X1 = [mesh(end,1), mesh(end,2)];
dy = X1(2)-X0(2);
dx = X1(1)-X0(1);

offsetAngle = -1*atan2(dy,dx);

%rotate everything by -theta to align to horizontal
% & offset everything to X0
[xC,yC,meshC] = rotAll2d(x,y,mesh, X0, offsetAngle);

% figure out wheterh avg mesh y is < 0. if so rotate 180 degrees
if mean([meshC(:,2);meshC(:,4)])<0
  [xC,yC,meshC] = rotAll2d(xC,yC,meshC, [0 0], pi);
  %update offset so it's zero again
  minX = min([meshC(:,1);meshC(:,3)]);
  minY = min([meshC(:,2);meshC(:,4)]);
  [xC,yC,meshC] = rotAll2d(xC,yC,meshC, [minX minY], 0);
end



%%tmp
%figure;
%subplot(2,1,1);
%hold all;
%plot(x,y,'r.');
%plot(mesh(:,1),mesh(:,2),'k');
%plot(mesh(:,3),mesh(:,4),'k');
%axis equal;
%
%subplot(2,1,2);
%hold all;
%plot(xC,yC,'r.');
%plot(meshC(:,1),meshC(:,2),'k');
%plot(meshC(:,3),meshC(:,4),'k');
%axis equal;
%pause

%-------------------------------------
function [xC,yC,meshC] = rotAll2d(x,y,mesh, X0, theta)
x= x-X0(1);
y= y-X0(2);
meshX1 = mesh(:,1)-X0(1);
meshY1 = mesh(:,2)-X0(2);
meshX2 = mesh(:,3)-X0(1);
meshY2 = mesh(:,4)-X0(2);


[xC,yC] = rot2d(x,y,theta);
xC = xC';
yC=yC';

[meshX1C, meshY1C]= rot2d(meshX1,meshY1,theta);
[meshX2C, meshY2C]= rot2d(meshX2,meshY2,theta);
meshC = [meshX1C',meshY1C',meshX2C',meshY2C'];

%-----------------------------
function [xC,yC] = rot2d(x,y,theta);
%xC, yC returned as 1xM vectors


R = [cos(theta) -sin(theta);...
     sin(theta) cos(theta)];

X = [x(:)';y(:)'];

Xc = R*X;

xC = Xc(1,:);
yC = Xc(2,:);





