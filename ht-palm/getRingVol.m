function [volRing nVoxRing] = getRingVol(vol,isoVal,xEdge,yEdge,zEdge)

xPix = xEdge(2) - xEdge(1);
yPix = yEdge(2) - yEdge(1);
zPix = zEdge(2) - zEdge(1);

vol = double(vol);

inRing = (vol>=isoVal);
nVoxRing = sum(inRing(:));
volRing = nVoxRing*xPix*yPix*zPix;

%-------------------------------------------------------
function [fv vol2 xEdge2 yEdge2 zEdge2] = closedRingVol(vol,isoVal,xEdge,yEdge,zEdge)

vol = double(vol);

%pad  the volume with zeros to close the isosurface 
volSz = size(vol);
vol2 = zeros(volSz+2);
vol2(2:end-1,2:end-1,2:end-1) = vol;
fv = isosurface(double(vol2),isoVal); 

%pad all edges  of the volume too
xPix = xEdge(2) - xEdge(1);
xEdge2 = [xEdge(1)-xPix, xEdge, xEdge(end)+xPix];
yPix = yEdge(2) - yEdge(1);
yEdge2 = [yEdge(1)-yPix, yEdge, yEdge(end)+yPix];
zPix = zEdge(2) - zEdge(1);
zEdge2 = [zEdge(1)-zPix, zEdge, zEdge(end)+zPix];

%---------------------------------------------------------
function [hP] = plotSurf(fv);

hP = patch(fv);
axis equal;
set(hP,'FaceColor','green','EdgeColor','none');

cameratoolbar;
camproj('perspective')
view(3);
hlight = camlight('headlight'); 
set(hP,'AmbientStrength',0.2,'SpecularStrength',0.2,'DiffuseStrength',.7);
set(hP,'BackFaceLighting','lit');
lighting gouraud
set(gcf,'Renderer','OpenGL')
camzoom(0.5)
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'Color','k');
set(gca,'Color','k');

axis(gca,'off');

