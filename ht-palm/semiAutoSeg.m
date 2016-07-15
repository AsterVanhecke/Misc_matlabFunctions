function [cellIm, recordSteps]= semiAutoSeg(im,param,h1)

if strcmp(param.edgemode,'valley')
  res = valleySegmentation2(im,param.threshminlevel,param.opennum,param.erodeNum,...
          param.edgeSigmaV,param.valleythresh1,param.valleythresh2,param.invertimage,param.thresFactorM);
elseif strcmp(param.edgemode,'sauvola')
  res = sauvola(im,param.localNeighbourhood,param.localThresh);
else
  error('Unrecognised edge mode.');
end
%remove isolated bg points
res=imfill(res,'holes');
%remove regions too small to be cells
res = bwareaopen(res,param.areaMin);

res0=res; %original image
resCur = res;%current image
resLast = res;% undo/ backup image

finished = false;
recordSteps=[];
if exist('h1','var')
  figure(h1);
else
  h1=figure;
end

fprintf('\n\n');
while ~finished
  %update fig
  figure(h1);
  labelIm = bwlabel(resCur);
  RGB  = label2rgb(labelIm,'jet','k','shuffle');
  imshow(RGB,'InitialMagnification',60);

  %options
  fprintf('Cell segmentation options:\n');
  fprintf(' (1) Split cells\n');
  fprintf(' (2) Join cells\n');
  fprintf(' (3) Delete cells\n');
  fprintf(' (4) Select good cells (deletes other cells)\n');
  fprintf(' (5) Smooth cells\n');
  fprintf(' (6) Delete edge cells\n');
  fprintf(' (7) Undo last change \n');
  fprintf(' (8) Discard all changes (reset to original image\n');
  fprintf(' (9) Finished segmentation\n');
  fprintf(' (10) Skip this FOV\n');
  %fprintf(' (10) Open cells\n');
  str = input('Choose 1-9: ','s');
  recordSteps=[recordSteps,str2num(str)];
  switch str
  case '1' %Split cells
    resLast = resCur;
    resCur = splitCell(resCur,param.areaMin,h1);
  case '2' % Join Cells
    resLast = resCur;
    resCur = joinCell(resCur);
  case '3'
    resLast = resCur;
    resCur = deleteCell(resCur);
  case '4'
    fprintf('\n1');
    resLast = resCur;
    resCur = selectGoodCell(resCur,h1);
  case '5' 
    resLast = resCur;
    fprintf('\nEnter radius of strel for smoothing via image closing and opening\n');
    fprintf('Default is %d, press Enter for default\n',param.smoothNum);
    str = input('Strel radius: ','s');
    if strcmp(str,'')
      smoothNum = param.smoothNum;
    else
      smoothNum = str2num(str);
    end
    resCur= smoothCell(resCur,smoothNum);
    recordSteps = [recordSteps,smoothNum];
  case '6'
    resLast = resCur;
    resCur = imclearborder(resCur);
  case '7'
    fprintf('Undid last change\n');
    resCur = resLast;
  case '8'
    fprintf('Undid all changes\n');
    resCur = res0;
  case '9'
    cellIm = resCur;
    finished = true;
  case '10' 
    resLast = resCur;
    resCur =resCur.*0;
    cellIm = resCur;
    finished = true;
  otherwise
    fprintf('\nInvalid choice\n\n');
  end
end

%----------------------------------------
function res = splitCell(res,minCellSize,h1)

figure(h1);
res = watershedSeedSeg4(res,minCellSize);
%----------------------------------------
function res0= joinCell(res)

[x,y] = ginputc_nice(size(res),2,'ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

%get the selected cells
labelIm = bwlabel(res);
nCell = numel(xI);
joinCell = zeros(nCell,1);
for ii = 1:nCell;
  labelIm = bwlabel(res);
  joinCell(ii) = labelIm(yI(ii),xI(ii));
end

joinCellIm= logical(zeros(size(res)));
joinCell = unique(joinCell);
for ii = 1:numel(joinCell)
  joinCellIm(find(labelIm==joinCell(ii)))=1;
end

%close the image until they are connected
isJoined=false;
ii=1;
while ~isJoined
  joinCellLabel = bwlabel(joinCellIm);
  if max(joinCellLabel(:))<=1
    isJoined=true;
  else
    joinCellIm = imclose(joinCellIm,strel('disk',ii));
    ii=ii+1;
  end
end

%copy the joined image back to the main image
res0=res;
res0(joinCellIm==1)=1;

%----------------------------------------
function goodCellIm= selectGoodCell(res,h1)
[x,y] = ginputc_nice(size(res),'ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

labelIm = bwlabel(res);
nCell = numel(xI);
goodCell = zeros(nCell,1);
for ii = 1:nCell;
  labelIm = bwlabel(res);
  goodCell(ii) = labelIm(yI(ii),xI(ii));
end

goodCellIm= logical(zeros(size(res)));
goodCell = unique(goodCell);
for ii = 1:numel(goodCell)
  goodCellIm(find(labelIm==goodCell(ii)))=1;
end

%%%----------------------------------------
%function closeCellIm= closeCell(res,closeNum)
%
%labelIm = bwlabel(res);
%zeroIm = logical(zeros(size(res)));
%oneCellIm = zeroIm;
%closeCellIm=zeroIm;
%
%nCell = max(labelIm(:));
%for ii = 1:nCell
%  oneCellIm = zeroIm;
%  oneCellIm(labelIm==ii)=1;
%  %smooth
%  oneCellIm= imclose(oneCellIm,strel('disk',closeNum));
%  closeCellIm(oneCellIm==1)=1;
%end
%%----------------------------------------
%function openCellIm= openCell(res,openNum)
%
%labelIm = bwlabel(res);
%zeroIm = logical(zeros(size(res)));
%oneCellIm = zeroIm;
%openCellIm=zeroIm;
%
%nCell = max(labelIm(:));
%for ii = 1:nCell
%  oneCellIm = zeroIm;
%  oneCellIm(labelIm==ii)=1;
%  %smooth
%  oneCellIm= imopen(oneCellIm,strel('disk',openNum));
%  openCellIm(oneCellIm==1)=1;
%end
%----------------------------------------
%----------------------------------------
function smoothCellIm= smoothCell(res,smoothNum)

labelIm = bwlabel(res);
zeroIm = logical(zeros(size(res)));
oneCellIm = zeroIm;
smoothCellIm=zeroIm;

nCell = max(labelIm(:));
for ii = 1:nCell
  oneCellIm = zeroIm;
  oneCellIm(labelIm==ii)=1;
  %smooth
  oneCellIm= imclose(oneCellIm,strel('disk',smoothNum));
  oneCellIm= imopen(oneCellIm,strel('disk',smoothNum));
  smoothCellIm(oneCellIm==1)=1;
end
%----------------------------------------
function badCellIm= deleteCell(res);
[x,y] = ginputc_nice(size(res),'ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

labelIm = bwlabel(res);
nCell = numel(xI);
badCell = zeros(nCell,1);
for ii = 1:nCell
  labelIm = bwlabel(res);
  badCell(ii) = labelIm(yI(ii),xI(ii));
end

badCellIm= res;
badCell = unique(badCell);
for ii = 1:numel(badCell)
  badCellIm(find(labelIm==badCell(ii)))=0;
end


