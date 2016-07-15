function [varargout] = deleteBadMesh2_Simple(srInMesh,limStruct, varargin)

nTrim = 10;
nGradSmooth = 20;

doOutputSrData = false;
doShowMesh = false;
doShowCDF = false;
doPlotTrend=false;
doFilterLengthTime = false;
doMaxTime = false;
nn = numel(varargin);
ii = 1;
while ii <= nn
   if strcmp(varargin{ii},'OutputSrData')
      doOutputSrData = true;
      ii = ii+1;
   elseif strcmp(varargin{ii},'ShowMesh')
      doShowMesh = true;
      htpalmSummaryFname= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'ShowCDF')
      doShowCDF=true;
      ii=ii+1;
   elseif strcmp(varargin{ii}, 'FilterLengthTime');
      doFilterLengthTime = true;
      tMax_filter = varargin{ii+1};
      nSigma = varargin{ii+2};
      ii = ii+3;
   elseif strcmp(varargin{ii},'MaxTime')
      doMaxTime = true;
      maxTime = varargin{ii+1};
      ii = ii+2;
   else 
      ii = ii+1;
   end
end

%lengthLimPc= 0.05;
%widthLimPc= 0.02;
%wiggleLimPc= 0.05;
%lengthPcForWiggle = 0.5;

lengthLimPc       =limStruct.lengthLimPc;
widthLimPc        =limStruct.widthLimPc;
wiggleLimPc       =limStruct.wiggleLimPc;
lengthPcForWiggle =limStruct.lengthPcForWiggle;
widthIntLimPc = limStruct.widthIntLimPc;

nCell = numel(srInMesh);
width  = zeros(nCell,1);
wiggle = width;
length = width;
wiggle = width;
tExpt = width;
widthInt = width;

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
   widthAll = sqrt( (mesh(:,1) - mesh(:,3)).^2 + (mesh(:,1) - mesh(:,3)).^2); 
   width(ii) =max(widthAll);
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
   widthInt(ii) = getInteralWidth(srInMesh{ii}.meshAligned);
end


lengthWiggleMin = percentileLim(length,lengthPcForWiggle);
shortCells = length<lengthWiggleMin;
wiggle(shortCells) = 0;

lengthLim(1)= percentileLim(length,lengthLimPc);
lengthLim(2)= percentileLim(length,1-lengthLimPc);
widthLim(1) = percentileLim(width,widthLimPc);
widthLim(2) = percentileLim(width,1-widthLimPc);
wiggleLim(1) = percentileLim(wiggle,wiggleLimPc);
wiggleLim(2) = percentileLim(wiggle,1-wiggleLimPc);
widthIntLim(1) = percentileLim(widthInt,widthIntLimPc);
widthIntLim(2) = percentileLim(widthInt,1-widthIntLimPc);


if doShowCDF
   figure;
   cellNo = 1:nCell;
   subplot(2,2,1);
   plot(cellNo,sort(length));
   hold all;
   xLowLim = find(sort(length) ==lengthLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(length) ==lengthLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,lengthLim(1),'rx');
   plot(xHighLim,lengthLim(2),'rx');
   ylabel('length');
   
   subplot(2,2,2);
   plot(cellNo,sort(width));
   hold all;
   xLowLim = find(sort(width) ==widthLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(width) ==widthLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,widthLim(1),'rx');
   plot(xHighLim,widthLim(2),'rx');
   ylabel('width');
   
   subplot(2,2,3);
   plot(cellNo,sort(wiggle));
   hold all;
   xLowLim = find(sort(wiggle) ==wiggleLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(wiggle) ==wiggleLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,wiggleLim(1),'rx');
   plot(xHighLim,wiggleLim(2),'rx');
   ylabel('wiggle');

   subplot(2,2,4);
   plot(cellNo,sort(widthInt));
   hold all;
   xLowLim = find(sort(widthInt) ==widthIntLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(widthInt) ==widthIntLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,widthIntLim(1),'rx');
   plot(xHighLim,widthIntLim(2),'rx');
   ylabel('widthInt');

end

if doShowMesh

   figure;
   movieData = load(htpalmSummaryFname);
   nMov = numel(movieData.fnameList);
   for ii = 1:nMov
      meshData = load(movieData.fnameList{ii}.phPreMeshName);
      srData= load(movieData.fnameList{ii}.srInMeshName);
      phIm = imread(movieData.fnameList{ii}.phPreImName);
      nCell = numel(meshData.cellList{1});
      
      hold off;
      imshow(phIm);
      hold all;
      for jj = 1:nCell 
         cellCurMesh = meshData.cellList{1}{jj};
         cellCurSr= srData.cellList{1}{jj};
         if ~isempty(cellCurMesh)&~isempty(cellCurSr)
	    
		 keyboard
	    lengthCur    = cellCurMesh.errorEstimate.length;
	    widthCur     = cellCurMesh.errorEstimate.width;
            wiggleCur    = cellCurMesh.errorEstimate.wiggleNorm;
	    widthIntCur = getInteralWidth(cellCurSr.meshAligned);

            if isGoodCell(cellCurSr,lengthCur,widthCur, wiggleCur, widthIntCur,lengthLim, widthLim, wiggleLim,widthIntLim)
               plotMesh(cellCurMesh,'g');
            else
               plotMesh(cellCurMesh,'r');
            end

         end
      end
      pause
   end
end


if doOutputSrData
   jj =1;
   for ii = 1:nCell
      tCur = tExpt(ii);
      cellCur = srInMesh{ii};
      if isGoodCell(cellCur,length(ii),width(ii), wiggle(ii), widthInt(ii),lengthLim, widthLim, wiggleLim,widthIntLim)
         srInMeshOut{jj} = cellCur;
         jj = jj+1;
      end
   end
   varargout = {srInMeshOut};
end

%------------------------------------------------------------------
function cellOk = isGoodCell(cellCur,length,width, wiggle, widthInt,lengthLim, widthLim, wiggleLim,widthIntLim)

badLength = (length<lengthLim(1) | length>lengthLim(2));
badWidth= (width<widthLim(1) | width>widthLim(2));
badWiggle = (wiggle<wiggleLim(1) | wiggle>wiggleLim(2));
badWidthInt=  (widthInt<widthIntLim(1) | widthInt>widthIntLim(2));
cellOk = ~badLength & ~badWidth &~badWiggle &~badWidthInt;


%------------------------------------------------------------------
function plotMesh(cellCur,varargin)

plotArg = varargin;

mesh = cellCur.mesh;
hold all;
plot(mesh(:,1),mesh(:,2),plotArg{:});
plot(mesh(:,3),mesh(:,4),plotArg{:});


%------------------------------------------------------------------
function maxWidthInternal = getInteralWidth(meshAligned)

widthInternal = abs(meshAligned(:,2)-meshAligned(:,4));
maxWidthInternal = max(widthInternal);

