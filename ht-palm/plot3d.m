function [hP fv] = plot3d(vol,isoVal)
fv = isosurface(double(vol),isoVal); 
hP = patch(fv);
axis equal;
set(hP,'FaceColor','green','EdgeColor','none');

cameratoolbar;
%camproj('perspective')
camproj('orthographic')
view(3);
hlight = camlight('headlight'); 
set(hP,'AmbientStrength',0.2,'SpecularStrength',0.2,'DiffuseStrength',.7);
set(hP,'BackFaceLighting','lit');
lighting gouraud
set(gcf,'Renderer','OpenGL')
%camzoom(0.5)
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'Color','k');
set(gca,'Color','k');

axis(gca,'off');
