function [cumSurface,cumVolume] = compute3DShape(X,Y)
% Compute surface and volume of curve rotating around x axis
% The calculation is done cumulatively for each segment
% inputs:
%    X: array of x coordinates
%    Y: array of y coordinates with length(X) = length(Y)
%
% outputs:
%    cumSurface: Cummulative surface area at each segment
%    cumVolume:  Cummulative volume at each segment
%    Note that length(cumSurface) = length(X)-1 = length(cumVolume)
% Author: Felix Schaber

% First compute geometrical centroid, arclength and area of the segments
areaCentroidY = areaCentroid(Y);
areaCentroidY = areaCentroidY(2:end); % a single line has no area
curveCentroidY = curveCentroid(Y);
curveCentroidY = curveCentroidY(2:end); % a single point is not a line

% Compute the cummulative length at each segment
[~,seglen] = arclength(X, Y, 'spline');
cumLen = cumsum(seglen);

% Integrate using the trapezoid rule
% TODO: Try different integration techniques
cumArea = cumtrapz(X,Y);
cumArea = cumArea(2:end); %Area of first point is 0

% Then use Pappus's Centroid Theorem to calculate surface and volume
cumSurface = 2*pi*curveCentroidY .* cumLen;  
cumVolume = 2*pi*areaCentroidY .* cumArea;
end

% Usage example
%Xsquare = [0 1 2 3 4 5];
%Ysquare = [5 5 5 5 5 5];
%[cumSurface,cumVolume] = compute3DShape(Xsquare',Ysquare');

function centroidsY = curveCentroid(Y)
% Helper function to compute the geometrical centroid of a discrete curve
% Note that this function assumes evenly spaced X points
centroidsY = cumsum(Y) ./ (1:length(Y))';
end

function centroidsY = areaCentroid(Y)
% Helper function to compute the centroid of the discrete area between
% curve and x axis for each slab
% Note this function assumes evenly spaced X points of width 1

% Compute centroid and area of each rectangular slab
slab_centroid = Y/2;
slab_area = Y; %Spacing between points is 1

% Overall centroid is sum of the individual centroids weighted by area
centroidsY = cumsum(slab_centroid .* slab_area) ./ cumsum(slab_area);
end