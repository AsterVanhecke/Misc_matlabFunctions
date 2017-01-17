function autoScatterWTmut(X1,Y1,X2,Y2,Xlab,Ylab)
% function to plot a scatterplot of wild type and mutant data, calculate
% the correlation coefficient, run autofig and set the legend.
% made to plot A LOT of scatterplots of the SIM data.
% does: hold on, plot, grid on, axis equal, x/ylabel, legend including 
% Spearman correlation coefficient and p-value.
% Author: Aster Vanhecke
hold on
plot(X1,Y1,'o','MarkerSize', 3,'MarkerFaceColor',[0.30 0.74 0.93], 'Color',[0 0.45 0.74],'LineWidth',0.5)
plot(X2,Y2,'^','MarkerSize', 3,'MarkerFaceColor',[0.93 0.69 0.13], 'Color',[0.85 0.325 0.098],'LineWidth',0.5)

[corCof1, pval1]=nanCorr(X1,Y1,'type','Spearman');
[corCof2, pval2]=nanCorr(X2,Y2,'type','Spearman');

% legendText={['Wild Type, Spearman Rho: ' num2str(corCof1) ', p= ' num2str(pval1)]...
%     ,['FtsW**I*, Spearman Rho: ' num2str(corCof2) ', p= ' num2str(pval2)]};
autofig(Xlab,Ylab, ['WT, Spearman: ' num2str(corCof1,2) ', p= ' num2str(pval1,2)],['Mut, Spearman: ' num2str(corCof2,2) ', p= ' num2str(pval2,2)])
end