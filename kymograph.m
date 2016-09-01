function M=kymograph(ZintAll,demo)
% function to plot kymograph
% kymograph(ZintAll,demo)
% inputs:
%--------
% ZintAll: cell(matlab) array, dimenions: cells(bacteria)-by-frames. Each
% cell(matlab) contains the intensity profile over the cell(bacteria)
% length, as a 2-by-Length vector.  the first row contains the length
% coordinates. (for now only supports equally spaced lengths with the same
% interval for each cell). The second row contains the intensity value for
% the respective lengths.
% demo: picks which type of kymograph to plot: 1 for kymograph per cell, 2
% for overall demograph (ordered by length, discards time info), 3:
% demograph per timepoint


% find the overal longest cell.
lengthz=cellfun(@length,ZintAll);
maxLength=max(max(lengthz));
maxLength=2*ceil((maxLength+1)/2); % increase value of maxLength to next even integer to prevent bugs
[num_cells, num_frames]=size(ZintAll);

switch demo
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
            x=1:maxFrame;x=x.*5;
            y=-maxLength*30:60:maxLength*30;
            plotKymo(M,x,y)
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

