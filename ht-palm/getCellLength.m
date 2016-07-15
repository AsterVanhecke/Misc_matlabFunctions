function cellLength = getCellLength(srInMesh)

nCell = numel(srInMesh);
cellLength = zeros(nCell,1);
for ii = 1:nCell
   cellLength(ii) = srInMesh{ii}.length;
end

