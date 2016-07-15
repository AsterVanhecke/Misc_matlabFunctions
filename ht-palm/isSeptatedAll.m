function srInMesh=isSeptatedAll(srInMesh,sliceWidth,sliceStep,lenRange,plotOn)

nCell = numel(srInMesh);


for ii=1:nCell
  ii
  cell0=srInMesh{ii};
  srInMesh{ii}=isSeptated(cell0,sliceWidth,sliceStep,lenRange,plotOn);
end
