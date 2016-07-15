function [l, t, w, d, cellData]=pltMTS_diameter_withCurv(srInMesh,varargin)

manClass=[];
doBadCell=false;
doManualClassify = false;
doTMax=false;
doFilterList=false;
tMax=Inf;
nArg=numel(varargin);
ii = 1;
while ii <= nArg
  if strcmp(varargin{ii},'ManualClassification')
    doManualClassify=true;
    manClass= varargin{ii+1};
    ii=ii+2;
  elseif strcmp(varargin{ii},'TMax')
    doTMax = true;
    tMax= varargin{ii+1};
    ii=ii+2;
  elseif strcmp(varargin{ii},'FilterList')
    doFilterList= true;
    filterList= varargin{ii+1};
    ii=ii+2;
  else
    ii=ii+1;
  end
end
scrsz = get(0,'ScreenSize');
h = figure('Position',[(scrsz(3)-1280)/3 (scrsz(4)-720)/3 1280 720],'color','w');

nCell = numel(srInMesh);
kk = 1;
for ii =1:nCell

  tExpt(kk) = srInMesh{ii}.tExpt_min;
  maxD = srInMesh{ii}.diameter.maxD;
  minD = srInMesh{ii}.diameter.minD;
  
  % Anna: get bacteria length and curvature
  bactLengthVect(ii) = srInMesh{ii}.bactLength;
  
  %-----------------------------CURVATURE---------------------------------%
  curvature = srInMesh{ii}.curvature;
  deriv1 = srInMesh{ii}.deriv1;
  deriv2 = srInMesh{ii}.deriv2;
  %septumPos = minD(1,3);
  septumPos1 = find(max(curvature(:, 1)));
  septumPos2 = find(max(curvature(:, 2)));
  distSeptumFromTerminal(ii) = minD(1, 1);
  curvatureAtSeptum1 = curvature(septumPos1, 1);
  curvatureAtSeptum2 = curvature(septumPos2, 2);
  if curvatureAtSeptum1 > curvatureAtSeptum2
    curvatureAtSeptumAllCellMax(ii) = curvatureAtSeptum1;
    curvatureAtSeptumAllCellMin(ii) = curvatureAtSeptum2;
  else
    curvatureAtSeptumAllCellMax(ii) = curvatureAtSeptum2;
    curvatureAtSeptumAllCellMin(ii) = curvatureAtSeptum1;
  end

    edg = srInMesh{ii}.mesh;
    edgX1 = edg(:, 1);
    edgY1 = edg(:, 2);
    edgX2 = edg(:, 3);
    edgY2 = edg(:, 4);   
 
    pointTg = septumPos1;
    set(0,'CurrentFigure',h);
    plot(edgX1, edgY1, 'bo'),
    hold on
    R = 1/curvatureAtSeptum1(pointTg);
    if deriv2(pointTg) > 0
       if  deriv1(pointTg) > 0
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
        if deriv1(pointTg) > 0
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
  
  %--------------------------END--CURVATURE-------------------------------%
  plotMinD(kk) = minD(1,2);
  plotMaxD(kk) = min(maxD(:,2));
  W(kk) = plotMinD(kk)/plotMaxD(kk);
  fov(kk) =srInMesh{ii}.fovNo;
  kk=kk+1;
end

cellOk=ones(nCell,1);

if tMax 
  tOk = tExpt<tMax;
  cellOk = cellOk & tOk(:);
end
if doFilterList
  cellOk = cellOk & filterList;
end
if doManualClassify
  %first exclude the bad cell
  isBad=logical(zeros(size(cellOk)));
  isBad(manClass.badCell) = 1;
  cellOk = cellOk & ~isBad;
  %annotate the septated cells as 'bad cells' 
  %and plot separately
  isDiv = logical(zeros(size(cellOk)));
  isDiv(manClass.postDiv)=1;
  %HACK - list the cells including diameter
  cellData = [(1:nCell)',tExpt(:),plotMinD(:),W(:)];
  cellData(isDiv,3:4)=0;
  cellData(~cellOk,:)=[];

  %END HACK

  %separate classification keeping the manual cells
  cellOk = cellOk & ~isDiv;

  %annotate the manual cells
  tDiv = tExpt(isDiv);
  tDiv = tDiv(:);
  dDiv = zeros(size(tDiv));
  lDiv = bactLengthVect(isDiv); 
  dDiv2 = [plotMinD(isDiv)', plotMaxD(isDiv)', W(isDiv)'];
  %add the non-segmented post-div cells
  pdExt = manClass.postDivExtra;
  nFov = size(pdExt,1);
  for ii = 1:nFov
    nExtra= pdExt(ii,3);
    if nExtra>=0
      tExt = pdExt(ii,2)*ones(nExtra,1);
      dExt = zeros(nExtra,1);
      dExt2 = zeros(nExtra,3);
      tDiv = [tDiv;tExt];
      dDiv = [dDiv;dExt];
      dDiv2 = [dDiv2;dExt2];
    end
  end

end

t= [tExpt(cellOk)';tDiv(:)];
w= [W(cellOk)';dDiv(:)];
wPriv = W(cellOk)';
d = [plotMinD(cellOk)';dDiv(:)];
a = pi*d.^2/4;
l = [bactLengthVect(cellOk)'; lDiv(:)];
lPriv = bactLengthVect(cellOk)';
tPriv = tExpt(cellOk);

lengthTemp = bactLengthVect(cellOk);
diamTemp = plotMinD(cellOk);
%-------------------------------------------------------------------------%
%  define rate (diameter costrinction rate respect to time and lenght)
%-------------------------------------------------------------------------%

[deriv1DiamVsTime, deriv2DiamVsTime] = dfdxc(tDiv,dDiv);
[deriv1DiamVsLeng, deriv2DiamVsLeng] = dfdxc(lPriv,dDiv);
%[deriv1DiamVsLeng, deriv2DiamVsLeng] = dfdxc(l,dDiv);

%--------------------------------------------------------------------------%


% Anna: plot septum diameter vs baceria length 
figure,
plot(bactLengthVect(cellOk), plotMinD(cellOk), 'm*')
title('Caulobacter septum diameter vs length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Min. diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

% Anna: plot curvature vs time
figure,
plot(tExpt(cellOk), curvatureAtSeptumAllCellMax(cellOk), 'm*')
title('Caulobacter septum max curvature vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Curvature []', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(plotMinD(cellOk), curvatureAtSeptumAllCellMax(cellOk), 'm*')
title('Caulobacter septum max curvature vs septum diameter', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Curvature []', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(distSeptumFromTerminal(cellOk), curvatureAtSeptumAllCellMax(cellOk), 'm*')
title('Caulobacter septum max curvature vs distance from the terminal point', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Curvature []', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(diff(lengthTemp), diamTemp(1: length(diamTemp)-1), 'b*')
title('Caulobacter septum diameter vs elongation', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Elongation [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(lengthTemp(1:length(lengthTemp)-1), diff(diamTemp), 'b*')
title('Variation of the septum diameter vs length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Diff Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(diff(lengthTemp), diff(diamTemp), 'b*')
title('Variation of the septum diameter respect to the length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Elongation [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Diff Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(tExpt(cellOk), lengthTemp, 'b*')
hold on
title('Caulobacter length vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

fitTimeLengthMod = polyfit(tExpt(cellOk), lengthTemp, 1);
%fitTimeLengthMod = griddedInterpolant(tExpt(cellOk), lengthTemp);
timeLegthMod = @(leng) leng./fitTimeLengthMod(1) - fitTimeLengthMod(2)./fitTimeLengthMod(1);
%timeFromLength = fitTimeLengthMod(lengthTemp);

plot(timeLegthMod(lengthTemp), lengthTemp, 'r--')
text(500, 3000, [ ' F(L)= ' num2str(round(fitTimeLengthMod(1))) '*t + ' num2str(round(fitTimeLengthMod(2))) ] , 'FontSize',12 )
text(500, 2500, '$$ F_{model}(l) = m*t + q $$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')
legend('data', 'fit')
hold off


%--------------------------------------------------------------------------%
figure;
subplot(2,1,1);
hold on;
plot(tExpt(cellOk),plotMinD(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv,'ko');
end
ylim([0 600]);

subplot(2,1,2);
hold on;
plot(tExpt(cellOk),W(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv,'ko');
end
ylim([0 1.2]);
title('Caulobacter septum waist ratio r/R vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
ylabel('Waist ratio r/R (septum diameter/max diameter)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );


figure;
subplot(2,1,1);
hold on;
plot(tExpt(cellOk),plotMinD(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv2(:,1),'ro');
end
ylim([0 600]);
legend('Caulobacters before division', 'Caulobacters post division')
title('Caulobacter septum diameter vs time - with post divisional cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Septum diameter [nm]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
subplot(2,1,2);
hold on;
plot(tExpt(cellOk),W(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv2(:,3),'ro');
end
ylim([0 1.2]);
legend('Caulobacters before division', 'Caulobacters post division')
title('Caulobacter septum waist ratio r/R vs time - with post divisional cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
ylabel('Waist ratio r/R (septum diameter/max diameter)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );


%try an area plot
t= [tExpt(cellOk)';tDiv(:)];
w= [W(cellOk)';dDiv(:)];
wPriv = W(cellOk)';
d = [plotMinD(cellOk)';dDiv(:)];
a = pi*d.^2/4;
l = [bactLengthVect(cellOk)'; lDiv(:)];
lPriv = bactLengthVect(cellOk)';

figure;
plot(t,a,'ko');
title('Caulobacter area vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Area [nm^2]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );

% Anna: fit area
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0 0],... % lower boundary paramters conditions
%                'Upper',[Inf Inf],... % upper boundary paramters conditions
%                'StartPoint',[1 1]);
% ft = fittype('a*(x)^n+b','problem','n','options',fo);
% 
% [fitRes, gof2] = fit(t, a, ft, 'problem', 2)
% 
% figure,
% plot( fitRes,'m')
% hold on
% plot(t, a, 'b*')
% 
% legend('Data','n=2')
% hold off

% Anna: adding fitting using seamus functions
[tc tg wMax fitOutTime]  = fitFeingold(t, w);
[tg wMax fitOutTime]  = fitFeingoldFixedStart(t, w, tc);
figure
plot(fitOutTime(:, 1), fitOutTime(:, 2), 'r')
hold on
plot(t, w, 'bo')
% C = 0.9;
% text(  4, 5, ' wMax.*sqrt(1- ((tDiv-tc)./(tg-tc)).^2)',...
%       'HorizontalAlignment','left',...
%       'VerticalAlignment','top',...
%       'BackgroundColor',[1 1 1],...
%       'FontSize',12 )
%line([350 0.5],[400 1],'Linestyle',':','color','k')
text(300, 1,[  num2str(wMax,'%.2f') ['$$\sqrt{1 - ( t - ' num2str(round(tc)) ' )/(' num2str(round(tg)) ' - '  num2str(round(tc)) '))^2}$$'] ],'interpreter','latex', 'FontSize',12 )
text(300, 1.2, '$$ F_{model}(t) = wMax*\sqrt{1 - ( ((t - t_c )/( t_g - t_c) )^2}$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Caulobacter septum waist ratio r/R vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with Feingold model', 'Caulobacter during division')
title('Caulobacter septum diameter vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off

% Anna: adding fitting using seamus functions
[tc tg wMax fitOut]  = fitFeingoldForDeriv( l(1: length(l)-1), diff(w), timeLegthMod);
[tg wMax fitOut]  = fitFeingoldFixedStartForDeriv( l(1: length(l)-1), diff(w), tc, timeLegthMod);
figure
plot(fitOut(:, 1), fitOut(:, 2), 'g')
hold on
time = timeLegthMod(l);
plot(time(1:length(time)-1), diff(w), 'bo')

text(350, -0.1,[  num2str(wMax,'%.2f') ['$$\sqrt{(1 - ( ( t - ' num2str(round(tc)) ' )/(' num2str(round(tg)) ' - '  num2str(round(tc)) ') )^2)^{-1}}*1/( ' num2str(round(tg)) ' - ' num2str(round(tc)) ')$$'] ],'interpreter','latex', 'FontSize',12 )
text(350, -0.3, '$$ F_{model}(t) = wMax*\sqrt{(1 -  ( (t - t_c )/( t_g - t_c) )^2 )^{-1}}*1/(t_g - t_c)$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Variation of septum diameters respect to the length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with Feingold model', 'Caulobacter during division')
title('Caulobacter septum diameter derivative vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum diameter/max diameter(along bacteria length) derivative','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off


% Anna: adding fitting using seamus functions
[tc tg wMax fitOutLength]  = fitFeingoldForLength( lPriv, wPriv, timeLegthMod);
[tg wMax fitOutLength]  = fitFeingoldFixedStartForLength( lPriv, wPriv, tc, timeLegthMod);
figure
plot(fitOutLength(:, 1), fitOutLength(:, 2), 'g')
hold on
time = timeLegthMod(lPriv);
plot(time, wPriv, 'bo')

text(300, 1,[  num2str(wMax,'%.2f') ['$$\sqrt{1 - ( t - ' num2str(round(tc)) ' )/(' num2str(round(tg)) ' - '  num2str(round(tc)) '))^2}$$'] ],'interpreter','latex', 'FontSize',12 )
text(300, 1.2, '$$ F_{model}(t) = wMax*\sqrt{1 - ( ((t - t_c )/( t_g - t_c) )^2}$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Septum diameters vs time (computed from length)', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with Feingold model', 'Caulobacter during division')
title('Caulobacter septum diameter derivative vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min] - Computed from length','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );

