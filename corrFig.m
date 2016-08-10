function [corCof]=corrFig(varargin)
% Function to examine correlation between vector X and Y.
% Plot scatterplot, draw a line, calculate and display some correlation
% coefficients. Maybe more.
% Getting input variable names doesnt work for expressions, e.g. Lg-Lc.
% As a temporary workaround you could define them as a single variable,
% e.g.: ElongAfterConstr=Lg-Lc; and use this as input for corrFig.

narginchk(2,3);
X=varargin{1};
Y=varargin{2};
if size(varargin,2)==3
    figTitle=varargin{3};
end
figure,
hold on
X=reshape(X, [], 1);
Y=reshape(Y, [], 1);
cfit=fit(X,Y,'poly1');
% plot(X,Y,'.')
plot(cfit,X,Y,'predfunc')
legend('Data', ['Fitted curve: Y= ' num2str(cfit.p1) '*X + ' num2str(cfit.p2)], 'prediction bounds')
corCof=corr(X,Y);
xlabel(inputname(1))
ylabel(inputname(2))
if size(varargin,2)==3
    title([figTitle ', correlation coefficient ' num2str(corCof)])
else
    title(['Correlation test between ' inputname(1) ' and ' inputname(2) ', correlation coefficient ' num2str(corCof)])
end

end