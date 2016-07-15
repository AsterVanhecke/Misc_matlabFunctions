function srInMeshOut = deleteBadMesh(srInMesh)

nTrim = 10;
nGradSmooth = 20;
nCell = numel(srInMesh);
lengthLimPc= 0.05;
widthLimPc= 0.02;
wiggleLimPc= 0.05;
lengthPcForWiggle = 0.5;


maxWidth = zeros(nCell,1);
wiggle = maxWidth;
length = maxWidth;

for ii = 1:nCell
   mesh = srInMesh{ii}.mesh;
   %calculate metrics:
   %  length
   %  width (max width probably)
   %  wiggliness (smoothed (probably exclude first and last) 
   %        use: sum(abs(dTheta)) ie wiggliness not total dTheta
   %        Only apply to long cells (>minLength), else zero
   %  

   %WIDTH
   width = sqrt( (mesh(:,1) - mesh(:,3)).^2 + (mesh(:,1) - mesh(:,3)).^2); 
   maxWidth(ii) =max(width);
   %LENGTH
   length(ii) = srInMesh{ii}.length;
   %WIGGLE
   centreLine = [mean([mesh(:,1),mesh(:,3)],2), mean([mesh(:,2),mesh(:,4)],2)];
   centreLineTrim = centreLine(nTrim+1:end-nTrim,:);
   nPts = numel(centreLineTrim);
   centreLineSmth = [smooth(centreLineTrim(:,1),nGradSmooth),smooth(centreLineTrim(:,2),nGradSmooth)];
   dx = diff(centreLineSmth(:,1));
   dy = diff(centreLineSmth(:,2));
   theta = atan2(dy,dx)*180/pi;
   dTheta = diff(theta);
   wiggle(ii) = sum(abs(dTheta)/nPts);
end

lengthWiggleMin = percentileLim(length,lengthPcForWiggle);
shortCells = length<lengthWiggleMin;
wiggle(shortCells) = 0;

length_lowLim = percentileLim(length,lengthLimPc);
length_highLim = percentileLim(length,1-lengthLimPc);
maxWidth_lowLim = percentileLim(maxWidth,widthLimPc);
maxWidth_highLim = percentileLim(maxWidth,1-widthLimPc);
wiggle_lowLim = percentileLim(wiggle,wiggleLimPc);
wiggle_highLim = percentileLim(wiggle,1-wiggleLimPc);

badLength = (length<length_lowLim | length>length_highLim);
badWidth= (maxWidth<maxWidth_lowLim | maxWidth>maxWidth_highLim);
badWiggle = (wiggle<wiggle_lowLim | wiggle>wiggle_highLim);

cellNo =1:nCell;

%figure;
%plot(cellNo,sort(length));
%hold all;
%xLowLim = find(sort(length) ==length_lowLim);
%xLowLim=xLowLim(1);
%xHighLim = find(sort(length) ==length_highLim);
%xHighLim=xHighLim(1);
%plot(xLowLim,length_lowLim,'rx');
%plot(xHighLim,length_highLim,'rx');
%ylabel('length');
%
%figure;
%plot(cellNo,sort(maxWidth));
%hold all;
%xLowLim = find(sort(maxWidth) ==maxWidth_lowLim);
%xLowLim=xLowLim(1);
%xHighLim = find(sort(maxWidth) ==maxWidth_highLim);
%xHighLim=xHighLim(1);
%plot(xLowLim,maxWidth_lowLim,'rx');
%plot(xHighLim,maxWidth_highLim,'rx');
%ylabel('maxWidth');
%
%figure;
%plot(cellNo,sort(wiggle));
%hold all;
%xLowLim = find(sort(wiggle) ==wiggle_lowLim);
%xLowLim=xLowLim(1);
%xHighLim = find(sort(wiggle) ==wiggle_highLim);
%xHighLim=xHighLim(1);
%plot(xLowLim,wiggle_lowLim,'rx');
%plot(xHighLim,wiggle_highLim,'rx');
%ylabel('wiggle');
%
%h1=figure
%for ii = 1:nCell
%   figure(h1);
%   mesh = srInMesh{ii}.mesh;
%   hold off;
%   plot(mesh(:,1),mesh(:,2),'k');
%   hold all;
%   plot(mesh(:,3),mesh(:,4),'k');
%   axis equal;
%
%   tString = ['Cell ',num2str(ii), ':'];
%   if badLength(ii)
%      tString = [tString, 'BAD LENGTH, '];
%   end
%   if badWidth(ii)
%      tString = [tString, 'BAD WIDTH, '];
%   end
%   if badWiggle(ii)
%      tString = [tString, 'BAD WIGGLE.'];
%   end
%   if  badLength(ii)||badWidth(ii)||badWiggle(ii)
%      tColor = 'r';
%   else
%      tColor = 'k';
%   end
%   title(tString,'Color',tColor);
%
%   length(ii)
%   wiggle(ii)
%   maxWidth(ii)
%
%   pause
%end

cellOk = ~badLength & ~badWidth &~badWiggle;
nCellOk = numel(find(cellOk));
srInMeshOut = cell(1,nCellOk);
jj =1;
for ii = 1:nCell
   if cellOk(ii)
      srInMeshOut{jj} = srInMesh{ii};
      jj = jj+1;
   end
end

