function M=kymograph(ZintAll,demogr,varargin)
% function to plot kymograph
% M=kymograph(ZintAll,demogr,varargin)
% inputs:
%--------
% ZintAll: cell(matlab) array, dimenions: cells(bacteria)-by-frames. Each
% cell(matlab) contains the intensity profile over the cell(bacteria)
% length, as a 2-by-Length vector.  the first row contains the length
% coordinates. (for now only supports equally spaced lengths with the same
% interval for each cell). The second row contains the intensity value for
% the respective lengths.
% demogr: picks which type of kymograph to plot: 1 for kymograph per cell, 2
% for overall demograph (ordered by length, discards time info), 3:
% length-ordered demograph per timepoint
% optional input:
% Tc: matrix of different times of interest for each cell(e.g. measured constriction
% onset) to be plotted on the kymographs (option demogr==1).
% [numberOfCells N]=size(Tc), with N the number of times you want to plot
% for each cell.
% T0: Integer. Time value to assign to first timepoint. either specify
% as 'T0' or 't0', or just insert a single value.
% Output: M: matrix for imageplot of last plot, most useful for the
% demograph.

narginchk(2,6)
inputi=1;

if demogr==1 % get optional input if needed
    while inputi <= length(varargin)
        if size(varargin{1,inputi},1)==size(ZintAll,1)
            Tc=varargin{1,inputi};
            inputi=inputi+1;
        elseif length(varargin{1,inputi})==2
            switch varargin{1,inputi}
                case {'tc','Tc'}
                    Tc=varargin{1,inputi+1};
                    inputi=inputi+2;
                case {'t0','T0'}
                    T0=varargin{1,inputi+1};
                    inputi=inputi+2;
                otherwise
                    error(['Input argument ' num2str(inputi+2) ' is invalid!'])
            end
        elseif isnumeric(varargin{1,inputi}) && length(varargin{1,inputi})==1
            T0=varargin{1,inputi};
            inputi=inputi+1;
        else
            error(['Input argument ' num2str(inputi+2) ' is invalid!'])
        end
    end
end

% find the overal longest cell.
lengthz=cellfun(@length,ZintAll);
maxLength=max(max(lengthz));
maxLength=2*ceil((maxLength+1)/2); % increase value of maxLength to next even integer to prevent bugs
[num_cells, num_frames]=size(ZintAll);

warning('off','MATLAB:legend:IgnoringExtraEntries') % supress legend warning

switch demogr
    case 1
        % Kymographs for each cell
        for cellIdx=1:num_cells
            maxFrame=find(lengthz(cellIdx,:),1,'last'); % find the last frame containing information
            M=NaN(maxLength,maxFrame);
            for frIdx=1:maxFrame
                if ~isempty(ZintAll{cellIdx,frIdx})
                    pos1=(maxLength/2)-floor(length(ZintAll{cellIdx,frIdx})/2);
                    M( pos1:(length(ZintAll{cellIdx,frIdx})+pos1-1) , frIdx ) = ...
                        ZintAll{cellIdx,frIdx}(2,:)';
                end
            end
            x=0:(maxFrame-1);x=x.*5;
            if exist('T0','var')
                x=x+T0;
            end
            y=-maxLength*30:60:maxLength*30;
            plotKymo(M,x,y)
            if exist('Tc','var')
                for ii=1:size(Tc,2)
                    plot([Tc(cellIdx,ii) Tc(cellIdx,ii)], [y(1) y(end)])
                    legend('midcell','TcF','TcX','Tz')
                end
            end
        end
    case 2
        % Demograph for all cells
        M=NaN(maxLength,num_cells*num_frames);
        for cellIdx=1:num_cells
            for frIdx=1:num_frames
                if ~isempty(ZintAll{cellIdx,frIdx})
                    pos1=(maxLength/2)-round(length(ZintAll{cellIdx,frIdx})/2);
                    M( pos1:(length(ZintAll{cellIdx,frIdx})+pos1-1) , (cellIdx-1)*num_frames+frIdx ) = ...
                        (ZintAll{cellIdx,frIdx}(2,:)./max(ZintAll{cellIdx,frIdx}(2,:)))';
                end
            end
        end
        M=lengthOrderKymo(M,lengthz);
        % Get rid of NaN crap
        firstIdx=find(~isnan(M(size(M,1)/2,:)),1);
        M=M(:,firstIdx:end);
        % TODO: normalize each intensity profile - maybe do this during loops
        % to make M
        y=-maxLength*30:60:maxLength*30;
        plotKymo(M,y)
    case 3
        % Demograph per timepoint
        for frIdx=1:num_frames
            M=NaN(maxLength,num_cells);
            for cellIdx=1:num_cells
                if ~isempty(ZintAll{cellIdx,frIdx})
                    pos1=(maxLength/2)-round(length(ZintAll{cellIdx,frIdx})/2);
                    M( pos1:(length(ZintAll{cellIdx,frIdx})+pos1-1) , cellIdx ) = ...
                        (ZintAll{cellIdx,frIdx}(2,:)/max(ZintAll{cellIdx,frIdx}(2,:)))';
                end
            end
            M=lengthOrderKymo(M,lengthz(:,frIdx));
            % Get rid of NaN crap
            firstIdx=find(~isnan(M(size(M,1)/2,:)),1);
            M=M(:,firstIdx:end);
            y=-maxLength*30:60:maxLength*30;
            plotKymo(M,y)
        end
        
    otherwise
        error('Invalid input. Should be 1, 2 or 3')
end

end
function intMatrix=lengthOrderKymo(M,lengthz)
% M is the matrix with lengths and corresponding intensities as columns for
% each cell
lengthsVec=reshape(lengthz',[1 numel(lengthz)]);
[~, orderIdx]= sort(lengthsVec);
intMatrix=M(:,orderIdx);
end

function plotKymo(varargin)
M=varargin{1,1};
M(1,1)=-0.05; % workaround for colormap: distinguish NaN from small values
if length(varargin)==3
    x=varargin{1,2};
    y=varargin{1,3};
    imW=0; imL=x(end);
    figure, imagesc(x,y,M);
elseif length(varargin)==2
    y=varargin{1,2};
    imW=0; imL=size(M,2);
    figure, imagesc(M, 'YData', [y(1) y(end)]);
else
    [imW, imL]=size(M);
    figure, imagesc(M);
end

try 
    load('\\files0\data\vanhecke\My Documents\MATLAB\cmap_parula_0black.mat');
    colormap(cmap)
catch
    colormap(parula)
end
hold on, plot([0 imL],[imW/2 imW/2],'r')
end

