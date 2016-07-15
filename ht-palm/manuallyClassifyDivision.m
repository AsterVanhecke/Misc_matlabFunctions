function manClassifyFile=manuallyClassifyDivision(srInMesh, htPalmSummary,tStart,blurRadius,srImLim)
%Classify cells (labelled and unlabelled) as pre-divisional, post-divisional, 
% or exclude (ambiguous/ bad segmentation)

manClassifyFile=[htPalmSummary(1:end-4),'.manualClassify.mat'];
%load each fov - plot the SR image
data = load(htPalmSummary);
fnameList = data.fnameList;

nCell = numel(srInMesh);
nFov = numel(fnameList);
pixSize = srInMesh{1}.pixSize;
manClass.postDiv=[];
manClass.badCell=[];
manClass.postDivExtra=[];;
h1 = figure;
%load(manClassifyFile,'manClass','ii');%IF YOU NEED TO RESTART FROM ii=nn uncomment this
for ii = 1:numel(fnameList)%and change 1 to nn
  ii
  srIm = makeSrIm(fnameList{ii}.srCorName,pixSize,blurRadius,srImLim);
  srImSat = imadjust(srIm,stretchlim(srIm),[0 1])*0.8;
  [cellFov cellFovNo] = getCellFov(srInMesh,ii);
  %plot the cells
  if ~isempty(cellFov{1})
    t = round(cellFov{1}{1}.tExpt_min)
  else
   t = getAcqTime(fnameList,ii,tStart)
  end

  [postDivC,postDivExtraC,badCellC]= selectCellState(srImSat,cellFov,pixSize,h1,cellFovNo);
  %save three lists
  % postDivisional [fovNo, srInMeshCellNo] LEFT CLICK
  % postDivisionalExtra [fovNo, time, nExtra] STILL LEFT CLICK
  % bad [fovNo, srInMeshCellNo] RIGHT CLICK
  [cellNoOut]=convertToHtpalmNo(cellFovNo,{postDivC,badCellC});
  manClass.postDiv = [manClass.postDiv;cellNoOut{1}(:)];
  manClass.badCell = [manClass.badCell;cellNoOut{2}(:)];
  manClass.postDivExtra = [manClass.postDivExtra; ii,t,postDivExtraC];

   %save after each FOV in case of crash
   save(manClassifyFile,'manClass','ii'); 
end

%------------------------------------------------------------
function [cellFov cellFovNo] = getCellFov(srInMesh,fovNo)

nCell = numel(srInMesh);
cellFov{1} ={};
cellFovNo=[];
kk=1;
for jj=1:nCell
  if srInMesh{jj}.fovNo ==fovNo
    cellFovNo(kk) = jj;
    cellFov{1}{kk} = srInMesh{jj};
    pixSize = srInMesh{jj}.pixSize;
    %modify the mesh so it plots!
    %cellFov{1}{kk}.mesh = cellFov{1}{kk}.mesh/pixSize;
    kk=kk+1;
  end
end

%------------------------------------------------------------
function  [postDiv,postDivExtra,badCell]= selectCellState(srImSat,cellFov,pixSize,h1,cellFovNo)

isOk = false;
while ~isOk
  %find srInMesh cells associated with the FOV, plot them
  h1 = dispcellsimple(cellFov,srImSat,'PixelSize',pixSize,'CellLabel',cellFovNo,'b-','LineWidth',1.5);
  % click the cells
  fprintf('Select cells: POST-left, BAD-right\n');
  [x,y,button] = ginputc_nice(size(srImSat),'FigHandle',h1,'ShowPoints',true,'ConnectPoints',false,'Color','b');
  hold all;
  plot(x(button==1),y(button==1),'o','MarkerEdgeColor', [0 0.6 0],'MarkerFaceColor', [0 0.8 0]);
  plot(x(button==3),y(button==3),'ro','MarkerFaceColor','r');
  str=input('Is the current selection ok? y/n (y): ','s');
  if strcmp(str,'y') || strcmp(str,'')
    isOk =true;
    %figure out the cell #s
    postDiv=[];
    postDivExtra=0;
    badCell=[];
    cellNo= findCellClick(cellFov,pixSize,x,y);%-1 means its not a selected cell
    for ii = 1:numel(cellNo)
      if button(ii)==1
        if cellNo(ii)==-1
          postDivExtra=postDivExtra+1;
        else
          postDiv= [postDiv,cellNo(ii)];
        end
      elseif button(ii)==3 && cellNo(ii)~=-1
        badCell= [badCell,cellNo(ii)];
      end
    end

  end
end
%-----------------------------------------------
function cellNo= findCellClick(cellFov,pixSize,x,y);%-1 means its not a selected cell

cellNo=[];
nCell=numel(cellFov{1});
nPt = numel(x);
for ii = 1:nPt
  xC=x(ii);
  yC=y(ii);
  isFound = false;
  jj=1;
  cellNo(ii)=-1;%ie not associated with any cell
  while ~isFound && jj<=nCell
    mesh = cellFov{1}{jj}.mesh;
    mesh = mesh/pixSize;
    meshPoly = [mesh(:,1:2);flipud(mesh(:,3:4))];%turn the mesh into a closed polygon
    inCell = inpolygon(xC,yC,meshPoly(:,1),meshPoly(:,2));
    %%DEBUG
    %hold off;
    %plot(meshPoly(:,1),meshPoly(:,2));
    %hold all;
    %plot(xC,yC,'x');
    %%DEBUG
    %pause
    if inCell
      cellNo(ii)=jj;
      isFound = true;
    else
      jj=jj+1;
    end
  end
end
%-----------------------------------------------------------------
function tAcq = getAcqTime(fnameList,fovNo,tStart)
data1 = loadjson(fnameList{1}.phPreImMeta);
tString1 = data1.FrameKey_0x2D_0_0x2D_0_0x2D_0.Time;
data = loadjson(fnameList{fovNo}.phPreImMeta);
tStringCur = data.FrameKey_0x2D_0_0x2D_0_0x2D_0.Time;

tAcq = getExptTime(tString1,tStringCur)+tStart;
%-------------------------------------------------------------------
function tExpt_min = getExptTime(tStringStart,tStringCur);
n0=datenum(tStringStart,'yyyy-mm-dd HH:MM:SS');
nC=datenum(tStringCur,'yyyy-mm-dd HH:MM:SS');

tExpt_min = (nC-n0)*1440;%days to minutes

%-------------------------------------------------------------------
function [cellListOut]=convertToHtpalmNo(cellFovNo,cellNoList);

  cellListOut ={};

  nList = numel(cellNoList);
  for ii = 1:nList
    lCur = cellNoList{ii};
    cellListOut{ii}= cellFovNo(lCur);
  end

%-------------------------------------------------------------------
function  srIm = makeSrIm(srCorName,pixSize,srBlurSigma,srImLim)
  
% load the SR analysis
load(srCorName,'srData');
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
x = srData.data(:,xcol);
y = srData.data(:,ycol);

srIm = stormHist2d(x,y,pixSize,'XLim',srImLim{1},'YLim',srImLim{2},'BlurSigma',srBlurSigma,'Normalize',0.0);
srIm = im2uint8(srIm);

%------------------------------------------------------------------------------
