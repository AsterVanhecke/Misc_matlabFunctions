function cellListOut = addExptTime(cellList,tExpt_min)
cellListOut = cellList;

for ii = 1:numel(cellListOut{1})
   if ~isempty(cellListOut{1}{ii})
      cellListOut{1}{ii}.tExpt_min = tExpt_min;
   end
end
