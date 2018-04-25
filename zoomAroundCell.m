function zoomAroundCell(cellData,cellIdx,pxlSize) % cellData = cellList.meshData{1,1}
if ~exist('pxlSize','var'),pxlSize=1;end
xlim([cellData{1,cellIdx}.box(1) , cellData{1,cellIdx}.box(1)+cellData{1,cellIdx}.box(3)]*pxlSize)
ylim([cellData{1,cellIdx}.box(2) , cellData{1,cellIdx}.box(2)+cellData{1,cellIdx}.box(4)]*pxlSize)
end