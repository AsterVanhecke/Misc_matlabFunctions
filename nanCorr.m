function [corCof, pval]=nanCorr(X,Y,varargin)
% Function to calculate the correlation coefficient, ignoring NaN-values.
X=reshape(X, [], 1);
Y=reshape(Y, [], 1);
% remove nan values from both X and Y. (order  is important)
Y=Y(~isnan(X));
X=X(~isnan(X));
X=X(~isnan(Y));
Y=Y(~isnan(Y));

switch length(varargin)
    case 0
        [corCof, pval]=corr(X,Y);
    case 2
        [corCof, pval]=corr(X,Y,varargin{1},varargin{2});
    case 4
        [corCof, pval]=corr(X,Y,varargin{1},varargin{2},varargin{3},varargin{4});
end
end