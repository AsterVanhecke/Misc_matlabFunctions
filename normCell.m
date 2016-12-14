function cellToNorm=normCell(cellToNorm)
% function to normalize diamFWHM and ZintAll
% Author: Aster Vanhecke
cellToNorm=cellfun(@normalize,cellToNorm,'UniformOutput',false,'ErrorHandler', @turnErrorIntoEmpty);
function x=normalize(x)
x(2,:)=x(2,:)/nanmean(x(2,:));
% actually this could be done by one line of hard-to-read code:
% ZintAllNormMean=cellfun(@(x) [x(1,:); x(2,:)/nanmean(x(2,:))],ZintAll,'UniformOutput',false,'ErrorHandler', @turnErrorIntoEmpty);