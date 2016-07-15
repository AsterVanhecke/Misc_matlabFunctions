function plotRingContraction(srInMeshAll,movieName,ringStatName,pixSize,ringRange,figName,tZero)

[stdR ,t,rCM, rCM_BW, stdR_BW , goodIm]= calculateRingRadius(srInMeshAll,movieName,ringStatName,pixSize);

ringCM = mean(ringRange,2);
poleRingDivide = -550;
isPolar = [ringCM<poleRingDivide]';
h1= figure;hold all;
plot(t(isPolar&goodIm),stdR(isPolar&goodIm),'.');
plot(t(~isPolar&goodIm),stdR(~isPolar&goodIm),'.');
plot(t(isPolar&goodIm),rCM(isPolar&goodIm),'.');
plot(t(~isPolar&goodIm),rCM(~isPolar&goodIm),'.');
plot(t(isPolar&goodIm),stdR_BW(isPolar&goodIm),'.');
plot(t(~isPolar&goodIm),stdR_BW(~isPolar&goodIm),'.');
plot(t(isPolar&goodIm),rCM_BW(isPolar&goodIm),'.');
plot(t(~isPolar&goodIm),rCM_BW(~isPolar&goodIm),'.');

legend('Std Polar','Std Ring','RCM polar', 'rCM ring','Std Polar BW','Std Ring BW','RCM polar BW', 'rCM ring BW');
xlabel('Time, min');
ylabel('Radial standard deviation, nm');

saveas(h1,'ring_radius_plot.fig');
save('ringRadiusResults.mat','t','isPolar','ringCM','ringRange','stdR','rCM', 'rCM_BW', 'stdR_BW' , 'goodIm');


% plot nice ring contraction info
load('ringRadiusResults.mat','t','isPolar','ringCM','ringRange','stdR','rCM', 'rCM_BW', 'stdR_BW' , 'goodIm');
tP = t(isPolar&goodIm);
tR = t(~isPolar&goodIm);
rP = stdR(isPolar&goodIm);
rR = stdR(~isPolar&goodIm);


tMin = 150;
tMax = Inf;
tOk = tR>=tMin & tR<=tMax;
tFit = tR(tOk);
rFit = rR(tOk);
[b,stats] = robustfit(tFit,rFit);
rFit_u = b(1) + b(2)*tFit;
b(2)
%pretty figure
hF=figure;hold all;
hax = gca;
hR = plot(tR+tZero,rR,'o');
hP = plot(tP+tZero,rP,'o');
hFit = plot(tFit+tZero, rFit_u ,'--');

set(hFit                          , ...
  'LineWidth'       , 3           ,...
  'Color', [204/255, 0, 0]       );
set(hR                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 5           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 .7 .7]  );
set(hP                         , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 5           , ...
  'MarkerEdgeColor' , 'none'      , ...
  'MarkerFaceColor' , [.75 .75 1] );


hLeg = legend('Mid-cell','Polar');
hXL = xlabel('Time  post-synchrony (min)');
hYL = ylabel('Radius (nm)');
title(['Grad= ',num2str(b(2))]);
%zlim([50,300]

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([ hXL, hYL], ...
    'FontName'   , 'AvantGarde');
set([hLeg, gca]             , ...
    'FontSize'   , 8           );
set([hXL, hYL]  , ...
    'FontSize'   , 10          );

set(gcf,...
      'Color','w',...
      'PaperPositionMode','auto');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );
saveas(gcf,figName);
