function [corCof]=corrFig(X,Y)
% Function to examine correlation between vector X and Y.
% Plot scatterplot, draw a line, calculate and display some correlation
% coefficients. Maybe more.
% Getting input variable names doesnt work for expressions, e.g. Lg-Lc.
% As a temporary workaround you could define them as a single variable,
% e.g.: ElongAfterConstr=Lg-Lc; and use this as input for corrFig.

figure,
hold on
X=reshape(X, [], 1);
Y=reshape(Y, [], 1);
cfit=fit(X,Y,'poly1');
% plot(X,Y,'.')
plot(cfit,X,Y,'predfunc')
title(['Correlation test between ' inputname(1) ' and ' inputname(2)])
legend('Data', ['Fitted curve: Y= ' num2str(cfit.p1) '*X + ' num2str(cfit.p2)], 'prediction bounds')
corCof=corr(X,Y);
xlabel(inputname(1))
ylabel(inputname(2))
title(['Correlation test between ' inputname(1) ' and ' inputname(2) ', correlation coefficient ' num2str(corCof)])

end