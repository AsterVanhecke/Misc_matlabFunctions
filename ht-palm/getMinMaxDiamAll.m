function [srInMeshOut] = getMinMaxDiamAll(srInMesh,showFit);

if ~exist('showFit','var')
  showFit = false;
else
  figure;
end

nCell = numel(srInMesh);

kk=1;
badCell = [];
for ii = 1:nCell
  ii
  try 
    [maxD, minD, c0] =getDiamMinMax(srInMesh{ii},showFit);
    srInMeshOut{kk}=srInMesh{ii};
    srInMeshOut{kk}.diameter.maxD = maxD;
    srInMeshOut{kk}.diameter.minD = minD;
    % Save the central position of the min diam
    srInMeshOut{kk}.c0 = c0;
    kk=kk+1;
    %xlim([0 3000]);
    %ylim([0 700]);
    %pause;
  catch ME
     fprintf(  '\n******* Cell skipped for following reason: *******\n');
     fprintf('%s\n',getReport(ME));
     fprintf(  '**************************************************\n');
     badCell = [badCell, ii];
  end
end
