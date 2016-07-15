function [corCof]=corrFig(X,Y)
% Function to examine correlation between vector X and Y.
% Plot scatterplot, draw a line, calculate and display some correlation
% coefficients. Maybe more.
figure,
hold on
X=reshape(X, [], 1);
Y=reshape(Y, [], 1);
cfit=fit(X,Y,'poly1');
% plot(X,Y,'.')
plot(cfit,X,Y,'predfunc')
title(['Correlation test between ' inputname(1) ' and ' inputname(2)])
corCof=corr(X,Y);
xlabel(inputname(1))
ylabel(inputname(2))
title(['Correlation test between ' inputname(1) ' and ' inputname(2) ', correlation coefficient ' num2str(corCof)])

end