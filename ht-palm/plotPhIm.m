function hF = plotPhIm(srInMesh,pixSize,cellNo,showOutline,scaleBarLen,satVal,varargin)

SBHEIGHT= 2;
EDGEOFFSET = 10;

doSortCells=false;
ii = 1;
while ii <= numel(varargin)
   if strcmp(varargin{ii},'SortCells')
      doSortCells=true;
      ii = ii + 1;
   else 
      ii=ii+1;
   end
end

nCell = numel(srInMesh);
cellLength = getCellLength(srInMesh);
if doSortCells
   [dum, order]=sort(cellLength(:,1));
else
   order = 1:nCell;
end

ii = order(cellNo)
cellData = srInMesh{ii};

phIm = cellData.phSubIm;
mesh = cellData.meshPix;
box = cellData.boxPix;
mesh= [mesh(:,1) - box(1), mesh(:,2) - box(2),mesh(:,3) - box(1),mesh(:,4) - box(2)];
mesh = mesh+1;

hold off;
hF = gcf;
phImSat = imadjust(phIm,stretchlim(phIm,[satVal, 1-satVal]),[]);

if scaleBarLen > 0
  lenPix = round(scaleBarLen/pixSize);
  sBarStart = [size(phImSat,1)- EDGEOFFSET - SBHEIGHT, size(phImSat,2)- EDGEOFFSET - lenPix]
  phImSat(sBarStart(1):sBarStart(1)+SBHEIGHT-1, sBarStart(2):sBarStart(2)+EDGEOFFSET-1) = 255;
end

imshow(phImSat);

if showOutline
  hold on;
  plot(mesh(:,1),mesh(:,2),'g-','LineWidth',2);
  plot(mesh(:,3),mesh(:,4),'g-','LineWidth',2);
end
