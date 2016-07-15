function [L,D]=projectToMeshMod(x,y,mesh,varargin)
% [L,D]=projectToMesh(X,Y,MESH,STEPLENGTH)
% This function projects each of the of points X,Y onto the centerline of
% the mesh MESH. For speed the function may optionally acquire STEPLENGTH
% variable, otherwise it is computed automatically. The function outputs
% the mesh coordinates L,D of the point.

if isempty(mesh) || isempty(x) || length(x)~=length(y), L=[]; D=[]; return; end
if isempty(varargin)
    stplng=0.5*sqrt(diff(mesh(:,1)+mesh(:,3)).^2+diff(mesh(:,2)+mesh(:,4)).^2); 
else
    stplng=varargin{1};
end
%calculate centreline
xCurve2 = (mesh(:,1)+mesh(:,3))/2;
yCurve2 = (mesh(:,2)+mesh(:,4))/2;
f = 5;
nPtSmooth = 3;
%*************************************
%This section of the code extends the first and last points "outside" the cell
% I dont see the point however - it just seems to cause errors? - misshapen cell outlines at ends
%
delta = f*median(stplng);
%Extend first and last points on centre line:
% compute the average vector of the last nPtSmooth (to avoid large changes in direction with very short step lengths at end
dxStart = xCurve2(1) - xCurve2(nPtSmooth);
dyStart = yCurve2(1) - yCurve2(nPtSmooth);
% normalise the length and set it equal to delta
dStart = delta * (dxStart^2+dyStart^2)^-.5 * [dxStart, dyStart];
xCurve2(1) = xCurve2(1) + dStart(1);
yCurve2(1) = yCurve2(1) + dStart(2);
%same for end
dxEnd= xCurve2(end) - xCurve2(end-nPtSmooth) ;
dyEnd= yCurve2(end) - yCurve2(end-nPtSmooth) ;
% normalise the length and set it equal to delta
dEnd = delta * (dxEnd^2+dyEnd^2)^-.5 * [dxEnd, dyEnd];
xCurve2(end) = xCurve2(end) + dEnd(1);
yCurve2(end) = yCurve2(end) + dEnd(2);
stplng = [-delta;cumsum(stplng)];
%**************************************
%stplng = [0;cumsum(stplng)];
%L is length along
L = zeros(size(x));
%D is distance from
D = zeros(size(x));
%I is segment #
I = zeros(size(x));
for i=1:length(x)
    %get perpendicular distance to each segment of curve
    %and the nearest point in each segment
    [dist,pn] = point2linedist(xCurve2,yCurve2,x(i),y(i));
    [tmp,j] = min(dist);
    L(i) = stplng(j)+distance(pn(j,:),[xCurve2(j) yCurve2(j)]);
    if (x(i)-xCurve2(j))*(yCurve2(j+1)-yCurve2(j))<(y(i)-yCurve2(j))*(xCurve2(j+1)-xCurve2(j))
        ori=1;
    else
        ori=-1;
    end
    D(i) = sqrt(dist(j))*ori;
    I(i) = j;
end
%%<DEBUG
%hold off;
%plot(xCurve2,yCurve2)
%hold all;
%plot(mesh(:,1),mesh(:,2))
%plot(mesh(:,3),mesh(:,4))
%plot(xCurve2(1:2),yCurve2(1:2),'k')
%plot(xCurve2(end-1:end),yCurve2(end-1:end),'r')
%axis equal
%pause
%hold off
%plot(L,D,'x')
%pause
%%DEBUG>

function c=distance(a,b)
c = (sum((a-b).^2))^0.5;
    
function [dist,pn] = point2linedist(xline,yline,xp,yp)
% point2linedist: distance,projections(line,point).
% A modification of SEGMENT_POINT_DIST_2D
% (http://people.scs.fsu.edu/~burkardt/m_src/geometry/segment_point_dist_2d.m)
% dist is the perpendicular distance
% Pn is the nearest point  on the line segment to Xp (perpendicularly)
dist = zeros(length(xline)-1,1);
pn = zeros(length(xline)-1,2);
p = [xp,yp];
for i=2:length(xline)
      p1 = [xline(i-1) yline(i-1)];
      p2 = [xline(i) yline(i)];
      if isequal(p1,p2)
          t = 0;
      else
          bot = sum((p2-p1).^2);
          t = (p-p1)*(p2-p1)'/bot;
          % if max(max(t))>1 || min(min(t))<0, dist=-1; return; end
          t = max(t,0);
          t = min(t,1);
      end
      pn(i-1,:) = p1 + t * ( p2 - p1 );
      dist(i-1) = sum ( ( pn(i-1,:) - p ).^2 );
end
