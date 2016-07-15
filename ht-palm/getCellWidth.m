function cell0 = getCellWidth(cell0,lStep,lWidth,histStep,showFit,nBreak);
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
diam_edge = zeros(nStep,1);
edgesLeft = zeros(nStep,2);
edgesRight = zeros(nStep,2);
fwhmLeft = zeros(nStep,2);
fwhmRight = zeros(nStep,2);
c0 = zeros(nStep, 2);
for ii = 1:nStep % Step along the bacteria length

  [slicePos{ii}, meshLim(ii,:), c0(ii, :), theta]=getCellSlice(x,y,mesh,cLineSmooth,len,l0(ii),lWidth,histStep);
  %TODO: find the FWHM (or blurred tophat fit)a -approach from the left, then the right
 [diam_fwhm(ii), diam_mesh(ii)]  = getCellFWHM(slicePos{ii},histStep,meshLim(ii,:));
 % [diam_edge(ii), edges, diam_fwhm(ii), diam_mesh(ii), xLeft, xRight,yLeft,yRight]  = getCellEdge(slicePos{ii},histStep,meshLim(ii,:));
  
  %Anna: rotate edges  
%   thetaRot = +theta;
%   R = [cos(thetaRot),-sin(thetaRot);...
%        sin(thetaRot),cos(thetaRot)];
%    %xEdgLeft = lMid(ii);
%    xEdgLeft = lWidth/2;
%    yEdgLeft = edges(1,1);
%    
%    %xEdgRight = lMid(ii);
%    xEdgRight = lWidth/2;
%    yEdgRight = edges(1,2);
%    
%   [XYrotLeft] = (R*[xEdgLeft;yEdgLeft])';
%   edgesLeft(ii,1) = XYrotLeft(:,1) + c0(1); %xEdgLeftRot
%   edgesLeft(ii,2) = XYrotLeft(:,2) + c0(2);  %yEdgLeftRot
%   
%   [XYrotRight] = (R*[xEdgRight;yEdgRight])';
%   edgesRight(ii,1) = XYrotRight(:,1) + c0(1);
%   edgesRight(ii,2) = XYrotRight(:,2) + c0(2);
%   
  %[XYrotLeftFWHM] = (R*[xLeft;yLeft])';
  %fwhmLeft(ii,1) = XYrotLeftFWHM(:,1) + c0(1); %xEdgLeftRot
%   [XYrotLeftFWHM] = (R*[xEdgLeft;xLeft])';
%   fwhmLeft(ii,1) = XYrotLeftFWHM(:,1) + c0(1);
%   fwhmLeft(ii,2) = XYrotLeftFWHM(:,2) + c0(2);  %yEdgLeftRot
  
%   [XYrotRightFWHM] = (R*[xEdgRight;xRight])';
%   %fwhmRight(ii,1) = XYrotRightFWHM(:,1) + c0(1);
%   fwhmRight(ii,1) = XYrotRightFWHM(:,1) + c0(1);
%   fwhmRight(ii,2) = XYrotRightFWHM(:,2) + c0(2);
  
  
  %%  IF its larger than the mesh width use the mesh width?
  %TODO - automatically find min (central) diameter and max (averaged somehow) width?-spline fit?
end

 idxPalmLocCrop=1;
 
 xL=mesh(:,1);
 yL=mesh(:,2);
 xR=mesh(:,3);
 yR=mesh(:,4);
 
% xLtemp=mesh(:,1);
% yLtemp=mesh(:,2);
% xRtemp=mesh(:,3);
% yRtemp=mesh(:,4);
% 
% %fitMeshLeft = polyfit(xLtemp, yLtemp, 3);
% fitMeshLeft=fit(xLtemp, yLtemp,'smoothingspline','SmoothingParam', 0.000001);
% % not more than 2nm precision expected
% xL = (min(xLtemp):4:max(xLtemp)-2);
% % yL = polyval(fitMeshLeft, xL);
% yL = fitMeshLeft(xL);
% 
% %fitMeshRight = polyfit(xRtemp, yRtemp, 3);
% fitMeshRight = fit(xRtemp, yRtemp, 'smoothingspline','SmoothingParam', 0.000001);
% % not more than 4nm precision expected
% xR = (min(xRtemp):4:max(xRtemp)-2);
% %yR = polyval(fitMeshRight, xR);
% yR = fitMeshRight(xR);
 
% for idxPalmLoc=1:length(x)
%     
%     idxYcorrPointL=find( abs(yL-y(idxPalmLoc))==min(abs(yL-y(idxPalmLoc))) );
%     idxYcorrPointR=find( abs(yR-y(idxPalmLoc))==min(abs(yR-y(idxPalmLoc))) );
%     
%     % if the point is inside the mesh line take it
%     if sum(x(idxPalmLoc)<xL(idxYcorrPointL))<1 && sum(x(idxPalmLoc)>xR(idxYcorrPointR))<1 && y(idxPalmLoc)<max([max(yL) max(yR)])&& y(idxPalmLoc)>min([min(yL) min(yR)])
%         xCrop(idxPalmLocCrop) = x(idxPalmLoc);
%         yCrop(idxPalmLocCrop) = y(idxPalmLoc);
%         idxPalmLocCrop = idxPalmLocCrop+1;
%     end
%    
% end
% 
%  if isempty(xCrop)
%      for idxPalmLoc=1:length(x)
% 
%     
%         idxYcorrPointL=find( abs(yL-y(idxPalmLoc))==min(abs(yL-y(idxPalmLoc))) );
%         idxYcorrPointR=find( abs(yR-y(idxPalmLoc))==min(abs(yR-y(idxPalmLoc))) );
%     
%         if sum(x(idxPalmLoc)>xL(idxYcorrPointL))<1 && sum(x(idxPalmLoc)<xR(idxYcorrPointR))<1 ... 
%                 && y(idxPalmLoc)<max([max(yL) max(yR)])&& y(idxPalmLoc)>min([min(yL) min(yR)])...
%             xCrop(idxPalmLocCrop) = x(idxPalmLoc);
%             yCrop(idxPalmLocCrop) = y(idxPalmLoc);
%             idxPalmLocCrop = idxPalmLocCrop+1;
%         end
%      end
%  end

 % DEBUG
% figure
% plot(x,y, 'r*');
% axis equal,
% hold on,
% plot(mesh(:,1), mesh(:,2),'y*' );
% plot(mesh(:,3), mesh(:,4),'y*' );
%  
% plot(xL, yL,'go' );
% plot(xR, yR,'go' );
% 
% plot(xCrop, yCrop, 'k*');
 
% 2D KDE of the PALM data
 % call the routine, which has been saved in the current directory 
 %   [bandwidth,density,X,Y]=kde2d([xCrop' yCrop'], 2^12);
  % plot the data and the density estimate
%     contour3(X,Y,density,50), hold on
%     plot(x,y,'r.','MarkerSize',5)
%idxPointsFilt = find( density>=mean(mean(density), 2) );

% ANNA "this is a Seamus length-line where each point is the center of the
%        bacteria slice" 
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
cell0.palmLoc = [x; y];
cell0.c0Pos = c0;
% cell0.diameter.diam_edge = diam_edge;
% cell0.edges = edges;

%  [C,h]=contour3(X,Y,density,3);

if showFit
    
  hold off;
  subplot(1,2,1)
  plot(lengthOut,diam_fwhm,'r-');
  hold all;
  plot(lengthOut,diam_mesh,'b-');
 % plot(lengthOut,diam_edge,'k-');
  %edgSmoothed = smooth(lengthOut, diam_edge, 'rloess', 4);
  
   %edgSmoothed2 = smooth(lengthOut, diam_edge, 'rloess', 4);
  %fitMod = griddedInterpolant(lengthOut, edgSmoothed, 'spline');
    % not more than 4nm precision expected
   % xx = (min(lengthOut):4:max(lengthOut));
   % edgSpline = fitMod(xx);
  %plot(lengthOut,edgSmoothed,'g-');
  %plot(xx,edgSpline,'c-');

  hold off;
 % legend('FWHM','Mesh', 'Edge', 'EdgeSmoothed','EdgeSpline');
 legend('FWHM','Mesh');
  xlabel('Length (nm)');
  ylabel('Diameter (nm)');
  ylim([0 700]);
  xlim([0 4000]);
  
%  
   subplot(1,2,2)
%   hold off;
  plot(x,y, 'k*');
  hold on,
%   plot(xCrop, yCrop, 'k*');
 
%   Ctemp = C(1,:);
%   idxCont = find(Ctemp<1);
%   idxDist = diff(idxCont);
%   longestContStartIdx=find(max(idxDist));
%   longestContStart = idxCont(longestContStartIdx);
%   longestContEnd = idxCont(longestContStartIdx)+idxDist(longestContStartIdx);
%   plot(C(1, longestContStart+1:longestContEnd-1), C(2, longestContStart+1:longestContEnd-1), 'c-');
  plot(mesh(:,1), mesh(:,2),'y-', 'LineWidth', 2);
  plot(mesh(:,3), mesh(:,4),'y-', 'LineWidth', 2);
  %plot(edgesLeft(:,1), edgesLeft(:,2), 'g-', 'LineWidth', 2);
 % plot(edgesRight(:,1), edgesRight(:,2), 'g-', 'LineWidth', 2);
%   plot(fwhmRight(:,1), fwhmRight(:,2), 'r-', 'LineWidth', 2);
%   plot(fwhmLeft(:,1), fwhmLeft(:,2), 'r-', 'LineWidth', 2); 
 % legend('PALM data','meshLeft', 'meshRight','edgeLeft','edgeRight', 'FWHMLeft','FWHMRight' )
 legend('PALM data','meshLeft', 'meshRight' )
  axis equal
  drawnow

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
%   xlabel('[nm]')
%   ylabel('[nm]')
%   hold all;
%   plot(xRotCrop,yRotCrop,'bx');
%   axis equal
%   subplot(2,1,2)
%   xCtr = meshLim(1)-histStep:histStep:meshLim(2)+histStep;
%   hist(slicePos,xCtr)
%    xlabel('Diameter[nm]')
%   ylabel('frequency (number of localization)')
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


function [diamEdge, edges, diamFWHM, diamMesh, posLeft, posRight,yLeft,yRight] = getCellEdge(slicePos,histStep,meshLim)
MINSTEP=1;%INTERPOLATION STEP - dont expect precision better than 1NM!
xCtr = meshLim(1)-histStep:histStep:meshLim(2)+histStep;
% create hist of the slice
nOut = hist(slicePos,xCtr);
diamMesh = meshLim(2)-meshLim(1);
[diamEdge, edges, diamFWHM, posLeft, posRight,yLeft,yRight]= getEdge(xCtr,nOut,meshLim,MINSTEP, slicePos,xCtr);


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
