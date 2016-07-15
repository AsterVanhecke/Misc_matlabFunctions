function [varargout] = deleteBadLength(srInMesh,tMax_filter,nSigma,varargin)

doOutputSrData = false;
doFilterLengthTime = false;
doFilterLengthTimeEqn = false;
doMaxTime = false;
doFilterLong = false;
doFilterShort=false;
doFilterVeryShort=false;
manClass=[];
doManualClassify = false;
nn = numel(varargin);
maxTime = Inf;
ii = 1;
while ii <= nn
   if strcmp(varargin{ii},'OutputSrData')
      doOutputSrData = true;
      ii = ii+1;
   elseif strcmp(varargin{ii},'ManualClassification')
    doManualClassify=true;
    manClass= varargin{ii+1};
    ii=ii+2;
   elseif strcmp(varargin{ii},'FilterLong')
      doFilterLong= true;
      ii = ii+1;
   elseif strcmp(varargin{ii},'FilterShort')
      doFilterShort= true;
      ii = ii+1;
   elseif strcmp(varargin{ii},'FilterVeryShort')
      doFilterVeryShort= true;
      ii = ii+1;
   elseif strcmp(varargin{ii}, 'FilterLengthTime');
      doFilterLengthTime = true;
      ii = ii+1;
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

nCell = numel(srInMesh);
width  = zeros(nCell,1);
len = width;
tExpt = width;

for ii = 1:nCell
   len(ii) = srInMesh{ii}.length;
   tExpt(ii) = srInMesh{ii}.tExpt_min;
end


if doMaxTime == true
  badTime = (tExpt>maxTime);
else
  badTime = zeros(nCell,1);
end

if doFilterLengthTimeEqn
 badLengthTimeEqn = len< (m*tExpt + c);
else
  badLengthTimeEqn= zeros(size(badTime));
end


cellOk = ~badTime & ~badLengthTimeEqn;

if doManualClassify
  %first exclude the bad cell
  isBad=logical(zeros(size(cellOk)));
  isBad(manClass.badCell) = 1;
  cellOk = cellOk & ~isBad;
  %annotate the septated cells as 'bad cells' 
  %and plot separately
  isDiv = logical(zeros(size(cellOk)));
  isDiv(manClass.postDiv)=1;
  cellOk = cellOk & ~isDiv;
  figure;
end

tGood = tExpt(cellOk);
lGood = len(cellOk);
t = tGood(tGood<tMax_filter);
l=  lGood(tGood<tMax_filter);
figure;

tCtr = 0:20:max(tGood);
lCtr = 1000:100:max(lGood);
ctrs={tCtr,lCtr}
figure;
n=hist3([tGood,lGood],'Ctrs',ctrs);
h = contourf(tCtr,lCtr,n');
colormap('gray');
cmap = colormap;
cmap = flipud(cmap);
colormap(cmap);
%set(h,'edgecolor','none');
hold all;



plot(tGood,lGood,'r.');
plot(t,l,'b.');

[b,stats] = robustfit(t,l)

std_robust = stats.robust_s;
lMax = [b(1)+std_robust*nSigma, b(2)];
lMin = [b(1)-std_robust*nSigma, b(2)];

plot(t,b(1)+b(2)*t,'k-');
plot(t,lMin(1)+lMin(2)*t,'k-');
plot(t,lMax(1)+lMax(2)*t,'k-');

if doFilterLengthTimeEqn
 plot(t,t*m + c,'r-');
end

%cells to keep
cellKeep= cellOk;
if doFilterLengthTime
  lMinOk = len>(lMin(1) + lMin(2)*tExpt);
  lMaxOk = len<(lMax(1) + lMax(2)*tExpt);
  cellKeep = lMinOk & lMaxOk & cellKeep;
end
if doFilterLong
  doFilterLong
  lOk = len> (b(1)+b(2)*tExpt);
  cellKeep = lOk&cellKeep;
end
if doFilterShort
  lOk = len< (b(1)+b(2)*tExpt);
  cellKeep = lOk& cellKeep;
end
if doFilterVeryShort
  lOk = len> (lMin(1) + lMin(2)*tExpt);
  cellKeep = lOk& cellKeep;
end

plot(tExpt(cellKeep),len(cellKeep),'g.');

plot([tMax_filter, tMax_filter],[min(l),max(l)],'k-');
xlabel('time (min)');
ylabel('length');

if doOutputSrData
   jj =1;
   nCell = numel(srInMesh)
   for ii = 1:nCell
      if cellKeep(ii)
         srInMeshOut{jj} = srInMesh{ii};
         jj = jj+1;
      end
   end
   varargout = {srInMeshOut,cellKeep};
end
%------------------------------------------------------------------
function cellOk = isGoodCell(cellCur,tCur,varargin)
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

len    = cellCur.len;

cellOk = true;;
if doFilterLengthTime
   lMin = lMinEq(1) + lMinEq(2)*tCur;
   lMax = lMaxEq(1) + lMaxEq(2)*tCur;
   badLen2 = len<lMin | len>lMax;
   cellOk = cellOk & ~badLen2;
end

if doFilterLengthTimeEqn
  lEqn = m*tCur+c;
  badLenEqn = len<lEqn;
  cellOk = cellOk & ~badLenEqn;
end

if doMaxTime
   badTime = (tCur>maxTime);
   cellOk = cellOk & ~badTime;
end

