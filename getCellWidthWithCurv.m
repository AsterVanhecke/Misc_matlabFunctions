function cell0 = getCellWidthWithCurv(cell0,lStep,lWidth,histStep,showFit,nBreak)
%Note - lStep is the distance along length for each slice
%       lWidth is the thickness of each slice - not the same
%
%%DEBUG
%global cLineSmooth;
%%DEBUG

if ~exist('showFit','var')
  showFit = false;
end
if ~exist('nBreak','var')
  nBreak= 10;
end

mesh = cell0.mesh;

cLine = zeros(size(mesh,1),2);
cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);
% smooth the centreLine
[x,y]=spline2d(cLine(:,1),cLine(:,2),nBreak);

cLineSmooth=[x(:),y(:)];

% Anna: compute bacteria length
bactLength = arclength(x, y, 'spline');
% Anna: compute curvature
[k1, deriv1Curv1, deriv2Curv1] = getCurvature(mesh(:,1), mesh(:,2));
[k2, deriv1Curv2, deriv2Curv2] = getCurvature(mesh(:,2), mesh(:,3));
k = zeros(length(k1), 2);
k(:,1) = k1;
k(:,2) = k2;
deriv1(:, 1)= deriv1Curv1;
deriv1(:, 2)=deriv1Curv2;
deriv2(:, 1)= deriv2Curv1;
deriv2(:, 2)=deriv2Curv2;

% Debug 
% disp('this is the distance computed with "arclength"')
% arclength(x, y, 'spline')
% disp('this is the distace computed with max min points')
% sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2)

%%DEBUG
%hold off
% plot(cLine(:,1),cLine(:,2),'r-');
% hold all;
% plot(x,y,'b-');
% plot(mesh(:,1),mesh(:,2),'g-');
% plot(mesh(:,3),mesh(:,4),'g-');
% axis equal
% legend('Raw','Smooth');
%pause;
%%DEBUG

len = sqrt((cLineSmooth(:,1)-cLineSmooth(1,1)).^2+(cLineSmooth(:,2)-cLineSmooth(1,2)).^2);% distances from the beginning of the mesh

%find points at which to crop the data and plot the histograms
l0 = 0:lStep:max(len);
lMid = l0(1:end-1)+diff(l0)/2;
nStep = numel(l0)-1;

xcol = findSRField(cell0.localizations.parInfo,'semantic','position in sample space in X'); 
ycol = findSRField(cell0.localizations.parInfo,'semantic','position in sample space in Y'); 
x = cell0.localizations.data(:,xcol);
y = cell0.localizations.data(:,ycol);
meshLim = zeros(nStep,2);
diam_fwhm = zeros(nStep,1);
diam_mesh = zeros(nStep,1);
for ii = 1:nStep
  [slicePos{ii},meshLim(ii,:)]=getCellSlice(x,y,mesh,cLineSmooth,len,l0(ii),lWidth,histStep);
  %TODO: find the FWHM (or blurred tophat fit)a -approach from the left, then the right
  [diam_fwhm(ii), diam_mesh(ii)]  = getCellFWHM(slicePos{ii},histStep,meshLim(ii,:));
  %%  IF its larger than the mesh width use the mesh width?
  %TODO - automatically find min (central) diameter and max (averaged somehow) width?-spline fit?
end
lengthOut = lMid';

cell0.diameter.slicePos = slicePos;
cell0.diameter.meshLim= meshLim;
cell0.diameter.lStep= lStep;
cell0.diameter.histStep= histStep;
cell0.diameter.length = lengthOut;
cell0.diameter.nBreak= nBreak;
cell0.diameter.diam_fwhm= diam_fwhm;
cell0.diameter.diam_mesh= diam_mesh;
% Anna 
cell0.bactLength = bactLength;
cell0.curvature = k;
cell0.deriv1 = deriv1;
cell0.deriv2 = deriv2;


if showFit

%   figure,
%   septumPos1 = find(max(curvature(:, 1)));
%   septumPos2 = find(max(curvature(:, 2)));
%   edg = srInMesh{ii}.mesh;
%    edgX1 = edg(:, 1);
%     edgY1 = edg(:, 2);
%     edgX2 = edg(:, 3);
%     edgY2 = edg(:, 4);   
%  
%     pointTg = septumPos1;
%     set(0,'CurrentFigure',h);
%     plot(edgX1, edgY1, 'bo'),
%     hold on
%     R = 1/curvatureAtSeptum1(pointTg);
%     if deriv2(pointTg) > 0
%        if  deriv1(pointTg) > 0
%             Xc = edgX1(pointTg) - R*sin( atan(dfdx(pointTg)) );
%             Yc = edgY1(pointTg) + R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%        else
%             Xc = edgX1(pointTg) - R*sin( atan(dfdx(pointTg)) );
%             Yc = edgY1(pointTg) + R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%        end
%         
%     else
%         if deriv1(pointTg) > 0
%             Xc = edgX1(pointTg) + R*sin( atan(dfdx(pointTg)) );
%             Yc = edgY1(pointTg) - R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%         else
%             Xc = edgX1(pointTg) + R*sin( atan(dfdx(pointTg)) );
%             Yc = edgY1(pointTg) - R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%             
%         end
%     end
    
  figure   
  hold off;
  plot(lengthOut,diam_fwhm,'r-');
  hold all;
  plot(lengthOut,diam_mesh,'b-');
  legend('FWHM','Mesh');
  xlabel('Length (nm)');
  ylabel('Diameter (nm)');
  ylim([0 700]);
  xlim([0 4000]);
end

%--------------------------------------------------------
function [slicePos,meshLim,c0,theta]=getCellSlice(x,y,mesh,cLine,len,l0,lWidth,histStep);
%find current mesh centre position, and +- lWidth/2
%probably best to smooth this!
cPt(:,1) = interp1(len,cLine(:,1),[l0, l0+lWidth/2, l0+lWidth]');
cPt(:,2) = interp1(len,cLine(:,2),[l0, l0+lWidth/2, l0+lWidth]');
c0 = cPt(2,:);%this is the midpoint

%get angle
dX = cPt(3,1)-cPt(1,1);
dY = cPt(3,2)-cPt(1,2);
theta = atan2(dY,dX);

%get meshMin Max distances
meshLim= getCellDist(mesh,c0,dY,dX);
if ~any(isnan(meshLim))
  %centre the xy data
  x = x-c0(1);
  y = y-c0(2);
  %rotate parallel to xAxis 
  thetaRot = -theta;
  R = [cos(thetaRot),-sin(thetaRot);...
       sin(thetaRot),cos(thetaRot)];
  [XYrot] = (R*[x';y'])';
  xRot = XYrot(:,1);
  yRot = XYrot(:,2);

  %crop
  isOk = xRot>=(-lWidth/2) & xRot<=(+lWidth/2) & yRot>=meshLim(1) & yRot<=meshLim(2);
  xRotCrop = xRot(isOk);
  yRotCrop = yRot(isOk);
  slicePos = yRotCrop;

  %DEBUG;
%   subplot(2,1,1)
%   hold off
%   plot(xRot,yRot,'rx');
%   hold all;
%   plot(xRotCrop,yRotCrop,'bx');
%   axis equal
%   subplot(2,1,2)
%   xCtr = meshLim(1)-histStep:histStep:meshLim(2)+histStep;
%   hist(slicePos,xCtr)
%   pause
  %%DEBUG;
else
    slicePos=[];
end



%---------------------------------------------------
function [meshLim] = getCellDist(mesh,c0,dY,dX)

%Get distance above
%rotate the gradient 90degerees to the cLine
perpLine_dX = -1*dY;
perpLine_dY = dX;
perpLine_m = perpLine_dY/perpLine_dX;
meshAbove= getDistLinePoly(mesh(:,1:2),c0,perpLine_m);
%Get distance below
meshBelow = getDistLinePoly(mesh(:,3:4),c0,perpLine_m);
meshLim = [-meshBelow, meshAbove];
%---------------------------------------------------
function distMax = getDistLinePoly(polyLine,c0,m);
%%DEBUG
%global cLineSmooth;
%%DEBUG
CELLRADMAX = 2000;
%Formula for (x1,y1) distance L along a straight line,
% given (x0,y0) and m 
y_along = @(y0,L,m,dirn) y0 +(dirn/abs(dirn))*m*L/sqrt(1+m^2);
x_along = @(x0,L,m,dirn) x0 +(dirn/abs(dirn))*L/sqrt(1+m^2);

%find intersection of line centred at c0(2,:) perpendicular to the centreline
%and the top mesh
%find the intersection
foundIsec=false;
distMax = [];
lineEnd = CELLRADMAX;
nMax = 2;
ii=1;
while foundIsec ==false
   X_line(1,:) = [x_along(c0(1),lineEnd,m,-1), y_along(c0(2),lineEnd,m,-1)];
   X_line(2,:) = [x_along(c0(1),lineEnd,m,+1), y_along(c0(2),lineEnd,m,+1)];
   [xIsec,yIsec] = intersections(polyLine(:,1),polyLine(:,2),X_line(:,1),X_line(:,2));
   %%DEBUG
   %hold off;
   %plot(cLineSmooth(:,1),cLineSmooth(:,2));
   %hold all;
   %plot(X_line(:,1),X_line(:,2));
   %plot(polyLine(:,1),polyLine(:,2));
   %if ~isempty(xIsec)
   %  plot(xIsec,yIsec,'ko');
   %end
   %axis equal
   %pause
   %%DEBUG

   if ~isempty(xIsec)
     distMax =sqrt((xIsec-c0(1)).^2+(yIsec-c0(2)).^2);
     distMax = min(distMax);
     foundIsec=true;
   else%increase the length of the search line (max cell radius) if no intersection found
     lineEnd = lineEnd*2;
     if ii ==nMax
       %Intersection between mesh and perpendicular to centreline not found!
       distMax=NaN;
       foundIsec=true;
     else
       ii=ii+1;
     end
   end
end

%-------------------------------------------------------------
function [xOut,yOut]=spline2d(x,y,nBreak);
nPt = numel(x);
nPtOut = 10*nPt;% upsample the output, but keep discrete

t=linspace(0,1,nPt);
tOut=linspace(0,1,nPtOut);
ppX = splinefit(t,x,nBreak);
ppY = splinefit(t,y,nBreak);
xOut = ppval(ppX,tOut);
yOut = ppval(ppY,tOut);


%--------------------------------------------
function [diamFWHM, diamMesh] = getCellFWHM(slicePos,histStep,meshLim)
MINSTEP=1;%INTERPOLATION STEP - dont expect precision better than 1NM!
xCtr = meshLim(1)-histStep:histStep:meshLim(2)+histStep;
nOut = hist(slicePos,xCtr);
[diamFWHM, posLeft, posRight,yLeft,yRight]= getFWHM(xCtr,nOut,meshLim,MINSTEP);
diamMesh = meshLim(2)-meshLim(1);

%%DEBUG
%hold off;
%bar(xCtr,nOut);
%hold all;
%plot(posLeft, yLeft,'ro');
%plot(posRight,yRight,'ro');
%pause
%%DEBUG
%%---------------------------------------------------
%function [diamMax diamMin]=findDiamMinMax(len,diam_fwhm);
%xSkipDist = 200;
%nBreak=10;
%diamMax=[];
%diamMin=[];
%
%%and chop off the last X nm becasue the ends are often imperfect
%isOk = len>xSkipDist& len<(len(end)-xSkipDist)
%len = len(isOk);
%diam_fwhm=diam_fwhm(isOk);
%
%%to find the max diameter - needs to be robust so use spline smoothing
%nPt = numel(len);
%ppDiam = splinefit(len,diam_fwhm,nBreak);
%diamSmth=ppval(ppDiam,len);
%
%[pks,idx]= findpeaks(diamSmth,'NPEAKS',2);
%if ~isempty(pks)
%  diamMax(1,:) = [len(idx(1)),pks(1)];
%  if numel(pks)>1
%    diamMax(2,:) = [len(idx(2)),pks(2)];
%  end
%end
%
%figure;
%hold all;
%plot(len,diamSmth,'k-')
%plot(len,diam_fwhm,'b-')
%plot(diamMax(:,1),diamMax(:,2),'ro');
%axis equal
%
%keyboard
%
