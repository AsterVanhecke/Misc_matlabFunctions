function [srInMeshOut,badCell]= getCellWidthAllWithCurv(srInMesh,lStep,lWidth,histStep,showFit)

if ~exist('showFit','var')
  showFit = false;
else
  figure
end

nCell = numel(srInMesh);

kk=1;
badCell = [];


% For each bacteria do this
for ii = 1:nCell
  ii
  try 
      % The output of getCellWidth is [cell0]
      % that is the diameter as a function of the length of the bacteria
      % and all the other info regarding each bacteria 
    srInMeshOut{kk} = getCellWidthWithCurv(srInMesh{ii},lStep,lWidth,histStep,showFit);
    kk=kk+1;
%      drawnow
  catch ME
     fprintf(  '\n******* Cell skipped for following reason: *******\n');
     fprintf('%s\n',getReport(ME));
     fprintf(  '**************************************************\n');
     badCell = [badCell, ii];
  end
end