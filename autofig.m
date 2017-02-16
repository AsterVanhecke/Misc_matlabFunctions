function autofig(Xlab,Ylab, varargin)
% function to automatically set figure lay-out, labels and (todo) limits.
% INPUT: Xlab, Ylab: abbreviation describing the 
% written for a jupyter notebook were A LOT of plots were made
% does: grid on, axis equal, x/ylabel, legend
% Author: Aster Vanhecke

grid on % I alwys like to have the grid as a reference

% set the legend
if isempty(varargin)
    legend('WT','Mut','Location', 'best')
else
    legend(varargin),
    legend('Location', 'best')
end

% write labels and get units
unitX=label(Xlab,'x');
unitY=label(Ylab,'y');

if strcmp(unitX,unitY) % if both axes have the same unit, ensure they have the same scale
   axis equal
end

end
function unit=label(lab,ax)
% function that evaluates INPUT "lab", and based on this sets the label of
% the axis (x or y, defined by INPUT "ax"). Then it determines the unit of
% this axis.
switch lab
    case 'dl'
        unit='um';
        if strcmp(ax,'x')
            xlabel('Elongation (µm)')
        else
            ylabel('Elongation (µm)')
        end
    case 'l'
        unit='um';
        if strcmp(ax,'x')
            xlabel('Length (µm)')
        else
            ylabel('Length (µm)')
        end
    case 't'
        unit='min';
        if strcmp(ax,'x')
            xlabel('Time (min)')
        else
            ylabel('Time (min)')
        end
    case 'dt'
        unit='min';
        if strcmp(ax,'x')
            xlabel('Duration (min)')
        else
            ylabel('Duration (min)')
        end
    case 'pd'
        unit='pd';
        if strcmp(ax,'x')
            xlabel('probability density')
        else
            ylabel('probability density')
        end
    case 'k'
        unit='k';
        if strcmp(ax,'x')
            xlabel('Exponent of elongation')
        else
            ylabel('Exponent of elongation')
        end
    case 'sav'
        unit='sav';
        if strcmp(ax,'x')
            xlabel('Surface area to volume ratio (µm^{-1})')
        else
            ylabel('Surface area to volume ratio (µm^{-1})')
        end
    case 'dsav'
        unit='sav';
        if strcmp(ax,'x')
            xlabel('Difference in surface area to volume ratio (µm^{-1})')
        else
            ylabel('Difference in surface area to volume ratio (µm^{-1})')
        end
end
end