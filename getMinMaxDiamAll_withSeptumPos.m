function [srInMeshOut] = getMinMaxDiamAll_withSeptumPos(srInMesh,showFit)

if ~exist('showFit','var')
  showFit = false;
else
  figure;
end

nCell = numel(srInMesh);

kk=1;
badCell = [];
scrsz = get(0,'ScreenSize');
h = figure('Position',[(scrsz(3)-1280)/3 (scrsz(4)-720)/3 1280 720],'color','w');
for ii = 1:nCell
  ii
  try 
    [maxD minD] =getDiamMinMax_getPosSeptum(srInMesh{ii},showFit);
    srInMeshOut{kk}=srInMesh{ii};
    srInMeshOut{kk}.diameter.maxD = maxD;
    srInMeshOut{kk}.diameter.minD = minD;
    kk=kk+1;
    %xlim([0 3000]);
    %ylim([0 700]);
    %pause;
    curvature = srInMesh{ii}.curvature;
    septumPos = minD(1,3);
  
    curvatureAtSeptum1 = curvature(septumPos, 1);
    curvatureAtSeptum2 = curvature(septumPos, 2);
    
    edg = srInMesh{ii}.mesh;
    edgX1 = edg(:, 1);
    edgY1 = edg(:, 2);
    edgX2 = edg(:, 3);
    edgY2 = edg(:, 4);
 
    pointTg = septumPos;
    set(0,'CurrentFigure',h);
    plot(edgX1, edgY1, 'bo'),
    hold on
    R = 1/curvatureAtSeptum1(pointTg);
    if d2fdx2(pointTg) > 0
       if  dfdx(pointTg) > 0
            Xc = edgX1(pointTg) - R*sin( atan(dfdx(pointTg)) );
            Yc = edgY1(pointTg) + R*cos( atan(dfdx(pointTg)) );
            plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
            axis equal
       else
            Xc = edgX1(pointTg) - R*sin( atan(dfdx(pointTg)) );
            Yc = edgY1(pointTg) + R*cos( atan(dfdx(pointTg)) );
            plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
            axis equal
       end
        
    else
        if dfdx(pointTg) > 0
            Xc = edgX1(pointTg) + R*sin( atan(dfdx(pointTg)) );
            Yc = edgY1(pointTg) - R*cos( atan(dfdx(pointTg)) );
            plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
            axis equal
        else
            Xc = edgX1(pointTg) + R*sin( atan(dfdx(pointTg)) );
            Yc = edgY1(pointTg) - R*cos( atan(dfdx(pointTg)) );
            plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
            axis equal
            
        end
    end
    drawnow
  catch ME
     fprintf(  '\n******* Cell skipped for following reason: *******\n');
     fprintf('%s\n',getReport(ME));
     fprintf(  '**************************************************\n');
     badCell = [badCell, ii];
  end
end
