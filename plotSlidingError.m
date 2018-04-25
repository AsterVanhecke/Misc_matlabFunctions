function varargout=plotSlidingError(x,y,smthWindow,varargin)
% Function to plot "sliding std" and/or "sliding sem"
% Calculates sliding average, std and SEM
% then uses the function "errorShade" to plot them.
% ErrorShade is adapted from Jean Yves Tinevez' msdanalyzer at:
% https://github.com/tinevez/msdanalyzer
% who adapted it from Rob Campbell code, at:
% http://www.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar/content/shadedErrorBar.m
% input: x,y,smthWindow, options
% x: x values
% y: y values 
% smthWindow
% options:
% 'std': add this string to plot standard deviation. It is possible to
% specify the color of the shaded area by adding a 1-b-3 vector with the
% RGB color code, e.g. [0.6 0.6 0.6], which is the default.
% 'sem': add this string to plot standard error of the mean. It is possible to
% specify the color of the shaded area by adding a 1-b-3 vector with the
% RGB color code, e.g. [0.1 0.1 0.1], which is the default.
% 'trim', N : remove (trim) a number of points (N, should be a discrete
% number) on the edge (std/sem is not meaningful at the edges, where only
% few points are used for the sliding average/std/sem.)

[x,idx]=sort(x);
y=y(idx);
x=reshape(x,[numel(x),1]);
y=reshape(y,[numel(y),1]);
plotSTD=false;
plotSEM=false;
colorsSEM=[0.1 0.1 0.1];
colorsSTD=[0.6 0.6 0.6];
trim=5;
%% parse input
argIdx=1;
while argIdx<=length(varargin)
    switch varargin{argIdx}
        case 'std'
            plotSTD=true;
            if length(varargin)>=argIdx+1
                if isnumeric(varargin{argIdx+1}) && numel(varargin{argIdx+1})==3 && length(varargin{argIdx+1})==3 %size(varargin{argIdx+1})==[1,3]
                    colorsSTD=varargin{argIdx+1};
                    argIdx=argIdx+1;
                end
            end
        case 'sem'
            plotSEM=true;
            if length(varargin)>=argIdx+1
                if isnumeric(varargin{argIdx+1}) && numel(varargin{argIdx+1})==3 && length(varargin{argIdx+1})==3
                    colorsSEM=varargin{argIdx+1};
                    argIdx=argIdx+1;
                end
            end
        case 'trim'
            trim=varargin{argIdx+1};
            argIdx=argIdx+2;
        otherwise
            error('Invalid input')
    end
    argIdx=argIdx+1;
end
%% Plotting
% get axes object
try
    ax=gca;
catch
    figure
    ax=gca;
end
% Create and plot sliding standard deviation
ySmth=smooth(y,smthWindow);
errorbars= sqrt(smooth((y-ySmth).^2,smthWindow));
errorbars(errorbars==0)=NaN;
if plotSTD
    H=errorShade(ax, x(1+trim:end-trim), ySmth(1+trim:end-trim), errorbars(1+trim:end-trim), colorsSTD, true);
end
% Create and plot sliding standard error of the mean
if plotSEM
    N=min([(1:length(errorbars))*2-1 ; 2*(length(errorbars):-1:1)-1 ; errorbars'.*0+smthWindow]);
    SEM_bars= errorbars./sqrt(N'); %
    errorShade(ax, x(1+trim:end-trim), ySmth(1+trim:end-trim), SEM_bars(1+trim:end-trim), colorsSEM, true);
end
if nargout>0
    varargout{1}=H;
end
end