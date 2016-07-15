function [varargout] = deleteBadMesh2(srInMesh,limStruct, varargin)

doOutputSrData = false;
doShowMesh = false;
doShowCDF = false;
doPlotTrend=false;
doFilterLengthTime = false;
doFilterLengthTimeEqn = false;
doMaxTime = false;
nn = numel(varargin);
maxTime = Inf;
tMax_filter=Inf;
nSigma=1.5;
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
   elseif strcmp(varargin{ii}, 'FilterLengthTimeEqn');
      doFilterLengthTimeEqn = true;
      m = varargin{ii+1};
      c = varargin{ii+2};
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
imForceLimPc         =limStruct.imForceLimPc;
nPixInOutLimPc       =limStruct.nPixInOutLimPc;
fitQLimPc          =limStruct.fitQLimPc;
widthIntLimPc = limStruct.widthIntLimPc;

nCell = numel(srInMesh);
width  = zeros(nCell,1);
wiggle = width;
length = width;
wiggle = width;
imForce = width;
nPixInOut = width;
fitQ = width;
tExpt = width;
widthInt = width;

for ii = 1:nCell
   %bug fixed 131217 correct length calculation
   length(ii) = getLen2(srInMesh{ii}.mesh);
   width(ii)     = srInMesh{ii}.errorEstimate.width;
   wiggle(ii)    = srInMesh{ii}.errorEstimate.wiggleNorm;
   imForce(ii)   = srInMesh{ii}.errorEstimate.imForceNorm;
   nPixInOut(ii) = srInMesh{ii}.errorEstimate.nPixInOutNorm;
   fitQ(ii)      = srInMesh{ii}.errorEstimate.fitquality;
   tExpt(ii) = srInMesh{ii}.tExpt_min;
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
imForceLim(1) = percentileLim(imForce,imForceLimPc);
imForceLim(2) = percentileLim(imForce,1-imForceLimPc);
nPixInOutLim(1) =percentileLim(nPixInOut,nPixInOutLimPc);   
nPixInOutLim(2) =percentileLim(nPixInOut,1-nPixInOutLimPc); 
fitQLim(1) = percentileLim(fitQ,fitQLimPc);
fitQLim(2) = percentileLim(fitQ,1-fitQLimPc);
widthIntLim(1) = percentileLim(widthInt,widthIntLimPc);
widthIntLim(2) = percentileLim(widthInt,1-widthIntLimPc);


if doShowCDF
   figure;
   cellNo = 1:nCell;
   subplot(3,2,1);
   plot(cellNo,sort(length));
   hold all;
   xLowLim = find(sort(length) ==lengthLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(length) ==lengthLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,lengthLim(1),'rx');
   plot(xHighLim,lengthLim(2),'rx');
   ylabel('length');
   
   subplot(3,2,2);
   plot(cellNo,sort(width));
   hold all;
   xLowLim = find(sort(width) ==widthLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(width) ==widthLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,widthLim(1),'rx');
   plot(xHighLim,widthLim(2),'rx');
   ylabel('width');
   
   subplot(3,2,3);
   plot(cellNo,sort(wiggle));
   hold all;
   xLowLim = find(sort(wiggle) ==wiggleLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(wiggle) ==wiggleLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,wiggleLim(1),'rx');
   plot(xHighLim,wiggleLim(2),'rx');
   ylabel('wiggle');

   subplot(3,2,4);
   plot(cellNo,sort(widthInt));
   hold all;
   xLowLim = find(sort(widthInt) ==widthIntLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(widthInt) ==widthIntLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,widthIntLim(1),'rx');
   plot(xHighLim,widthIntLim(2),'rx');
   ylabel('widthInt');

   subplot(3,2,5);
   plot(cellNo,sort(nPixInOut));
   hold all;
   xLowLim = find(sort(nPixInOut) ==nPixInOutLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(nPixInOut) ==nPixInOutLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,nPixInOutLim(1),'rx');
   plot(xHighLim,nPixInOutLim(2),'rx');
   ylabel('nPixInOut');

   subplot(3,2,6);
   plot(cellNo,sort(fitQ));
   hold all;
   xLowLim = find(sort(fitQ) ==fitQLim(1));
   xLowLim=xLowLim(1);
   xHighLim = find(sort(fitQ) ==fitQLim(2));
   xHighLim=xHighLim(1);
   plot(xLowLim,fitQLim(1),'rx');
   plot(xHighLim,fitQLim(2),'rx');
   ylabel('fitQ');
end



if doFilterLengthTime || doFilterLengthTimeEqn
   badLength = (length<lengthLim(1) | length>lengthLim(2));
   badWidth= (width<widthLim(1) | width>widthLim(2));
   badWiggle = (wiggle<wiggleLim(1) | wiggle>wiggleLim(2));
   badImForce =  (imForce<imForceLim(1) | imForce>imForceLim(2));
   badNPixInOut =  (nPixInOut<nPixInOutLim(1) | nPixInOut>nPixInOutLim(2));
   badFitQ =  (fitQ<fitQLim(1) | fitQ>fitQLim(2));
   badWidthInt=  (widthInt<widthIntLim(1) | widthInt>widthIntLim(2));
   if doMaxTime == true
      badTime = (tExpt>maxTime);
   else
      badTime = zeros(size(badLength));
   end

   if doFilterLengthTimeEqn
     badLengthTimeEqn = length< (m*tExpt + c);
   else
      badLengthTimeEqn= zeros(size(badLength));
   end

  
   cellOk = ~badLength & ~badWidth &~badWiggle & ~badImForce & ~badNPixInOut & ~badFitQ & ~badWidthInt & ~badTime & ~badLengthTimeEqn;

   t = tExpt(cellOk);
   l = length(cellOk);
   t = t(t<tMax_filter);
   l=  l(t<tMax_filter);
   figure;
   hold all;
   t = tExpt(cellOk);
   l = length(cellOk);
   plot(tExpt,length,'r.');
   plot(t,l,'b.');

   if doFilterLengthTime
     [b,stats] = robustfit(t,l)

     std_robust = stats.robust_s;
     lMax = [b(1)+std_robust*nSigma, b(2)];
     lMin = [b(1)-std_robust*nSigma, b(2)];

     plot(t,b(1)+b(2)*t,'k-');
     plot(t,lMin(1)+lMin(2)*t,'k-');
     plot(t,lMax(1)+lMax(2)*t,'k-');
   end

   if doFilterLengthTimeEqn
     plot(t,t*m + c,'r-');
   end
    
   cellLength = length;

   plot([tMax_filter, tMax_filter],[min(l),max(l)],'k-');
   xlabel('time (min)');
   ylabel('length');
end


if doShowMesh
   kk =1;
   goodCellVarArg = {};
   if doFilterLengthTime
      goodCellVarArg = {goodCellVarArg{:},lMin,lMax};
   end

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
            tCur = tExpt(kk)
            kk=kk+1;
            if isGoodCell(cellCurSr,lengthLim, widthLim, wiggleLim,widthIntLim,tCur,goodCellVarArg{:})
               plotMesh(cellCurMesh,'g');
            else
               plotMesh(cellCurMesh,'r');
            end

         end
      end
      pause
   end
end

goodCellList=zeros(nCell,1);
if doOutputSrData
   goodCellVarArg = {};
   if doFilterLengthTime
      goodCellVarArg = {goodCellVarArg{:},'LengthTime',lMin,lMax};
   end
   if doFilterLengthTimeEqn
      goodCellVarArg = {goodCellVarArg{:},'LengthTimeEqn',m,c};
   end
   if doMaxTime
      goodCellVarArg = {goodCellVarArg{:},'MaxTime',maxTime};
   end
   jj =1;
   nCell = numel(srInMesh);
   for ii = 1:nCell
      tCur = tExpt(ii);
      cellCur = srInMesh{ii};
      if isGoodCell(cellCur,lengthLim, widthLim, wiggleLim,widthIntLim,tCur,goodCellVarArg{:})
         srInMeshOut{jj} = cellCur;
         goodCellList(ii) = 1;
         jj = jj+1;
      end
   end
   varargout = {srInMeshOut,goodCellList};
end
%------------------------------------------------------------------
function cellOk = isGoodCell(cellCur,lengthLim, widthLim, wiggleLim,widthIntLim,tCur,varargin)
doFilterLengthTime = false;
doFilterLengthTimeEqn = false;
doMaxTime = false;
ii = 1;
while ii <= numel(varargin)
   if strcmp(varargin{ii},'LengthTime')
      doFilterLengthTime = true;
      lMinEq = varargin{ii+1};
      lMaxEq = varargin{ii+2};
      ii = ii +3;
   elseif strcmp(varargin{ii},'LengthTimeEqn')
      doFilterLengthTimeEqn = true;
      m = varargin{ii+1};
      c = varargin{ii+2};
      ii = ii +3;
   elseif strcmp(varargin{ii},'MaxTime')
      doMaxTime = true;
      maxTime = varargin{ii+1};
      ii = ii +2;
   else 
      ii = ii+1;
   end
end

length    = getLen2(cellCur.mesh);
width     = cellCur.errorEstimate.width;
wiggle    = cellCur.errorEstimate.wiggleNorm;
widthInt = getInteralWidth(cellCur.meshAligned);

badLength = (length<lengthLim(1) | length>lengthLim(2));
badWidth= (width<widthLim(1) | width>widthLim(2));
badWiggle = (wiggle<wiggleLim(1) | wiggle>wiggleLim(2));
badWidthInt=  (widthInt<widthIntLim(1) | widthInt>widthIntLim(2));
cellOk = ~badLength & ~badWidth &~badWiggle &~badWidthInt;

if doFilterLengthTime
   lMin = lMinEq(1) + lMinEq(2)*tCur;
   lMax = lMaxEq(1) + lMaxEq(2)*tCur;
   badLen2 = length<lMin | length>lMax;
   cellOk = cellOk & ~badLen2;
end

if doFilterLengthTimeEqn
  lEqn = m*tCur+c;
  badLenEqn = length<lEqn;
  cellOk = cellOk & ~badLenEqn;
end

if doMaxTime
   badTime = (tCur>maxTime);
   cellOk = cellOk & ~badTime;
end


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
%-------------------------------------------------------------
function len = getLen2(mesh)

centreLine = [mean([mesh(:,1),mesh(:,3)],2), mean([mesh(:,2),mesh(:,4)],2)];

dx = diff(centreLine(:,1));
dy = diff(centreLine(:,2));

%LENGTH
len = sum(sqrt((dx.^2+dy.^2)));
