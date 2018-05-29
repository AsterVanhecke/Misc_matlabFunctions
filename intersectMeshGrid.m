function [xIntersect, yIntersect]=intersectMeshGrid(mesh,x,y)
% function to asign L and D cell coordinates to x and y pixel grids.
% input:
% mesh is N-by-4. Each row of mesh holds [x1 y1 x2 y2] defining the segment
% border. 
% x an y are row vectors (1-by-M), where each column of x and y holds a value for
% the pixel grid.
% Output:
% xIntersect and yIntersect are N-by-M matrices,

% Convert mesh to slope m and offset q of lines defining segment borders.
% m = (y2-y1)/(x2-x1);
m = (mesh(:,4)-mesh(:,2))./(mesh(:,3)-mesh(:,1));
% q = y1 - m*x1
q = mesh(:,2) - m.*mesh(:,1);

% Each row corresponds to a segment.
% Each column corresponds to a pixel grid line.

% Find intersection between vertical pixel grid and segment borders:
% x = a
% y = m*x + q
% solution: y = m*a + q
yIntersect = m.*x+q;

% Find intersection between horizontal pixel grid and segment borders:
% y = b
% y = m*x + q
% solution: x = (b-q)/m
xIntersect = (y-q)./m;

% points of intersection are:
% x,yIntersect
% xIntersect,y
