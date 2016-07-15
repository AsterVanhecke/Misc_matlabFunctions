%------------------------------------------------------
function [srInMesh htPalmConfigData]= loadAllData(htPalmSummary,boxExtraSpace)

if ~exist('boxExtraSpace','var')
   boxExtraSpace  = 0;
end

srInMesh = {};
data = load(htPalmSummary);
pixNmPh =data.configData.pixUmPh*1e3;

fnameList = data.fnameList;
for ii = 1:numel(fnameList)
   cellListTmp1 = load(fnameList{ii}.srInMeshName);
   cellListTmp2 = cellListTmp1.cellList{1};
   for jj = 1:numel(cellListTmp2)
      if ~isempty(cellListTmp2{jj})
        cellListTmp2{jj}.fovNo = ii;
        phName  = (fnameList{ii}.phPreImName);
        [phSubIm,boxPix,meshPix]= getSubIm(cellListTmp2{jj},phName, pixNmPh,boxExtraSpace);

        cellListTmp2{jj}.phSubIm = phSubIm;
        cellListTmp2{jj}.boxPix = boxPix;
        cellListTmp2{jj}.meshPix = meshPix;
        cellListTmp2{jj}.tExpt_min= fnameList{ii}.tExpt_min;

        srInMesh = {srInMesh{:},cellListTmp2{jj}};
      end
   end
end

htPalmConfigData = data.configData;

%----------------------------------------------
function [phSubIm,boxPix,meshPix]= getSubIm(cellData,movieName,pixNmPh,boxExtraSpace)
phIm = imread(movieName);
boxPix = cellData.box/pixNmPh;
meshPix = cellData.mesh/pixNmPh;

boxPix(1) = boxPix(1)-boxExtraSpace;
boxPix(2) = boxPix(2)-boxExtraSpace;
boxPix(3) = boxPix(3)+2*boxExtraSpace;
boxPix(4) = boxPix(4)+2*boxExtraSpace;


phSubIm = phIm(boxPix(2):boxPix(2)+boxPix(4)-1,boxPix(1):boxPix(1)+boxPix(3)-1);

