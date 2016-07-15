function [z, rad,initguess] = fitCircleParametric(x, errorRatioYX, varargin)
% Fit a circle with asymettric errors in X & Y 
% rescale the data and fit an ellipse
%FITCIRCLE    least squares circle fit
%   
%   [Z, R] = FITCIRCLE(X) fits a circle to the N points in X minimising
%   geometric error (sum of squared distances from the points to the fitted
%   circle) using nonlinear least squares (Gauss Newton)
%       Input
%           X : 2xN array of N 2D points, with N >= 3
%       Output
%           Z : center of the fitted circle
%           R : radius of the fitted circle
%
%   [Z, R] = FITCIRCLE(X, 'linear') fits a circle using linear least
%   squares minimising the algebraic error (residual from fitting system
%   of the form ax'x + b'x + c = 0)
%
%   [Z, R] = FITCIRCLE(X, Property, Value, ...) allows parameters to be
%   passed to the internal Gauss Newton method. Property names can be
%   supplied as any unambiguous contraction of the property name and are 
%   case insensitive, e.g. FITCIRCLE(X, 't', 1e-4) is equivalent to
%   FITCIRCLE(X, 'tol', 1e-4). Valid properties are:
%
%       Property:                 Value:
%       --------------------------------
%       maxits                    positive integer, default 100
%           Sets the maximum number of iterations of the Gauss Newton
%           method
%
%       tol                       positive constant, default 1e-5
%           Gauss Newton converges when the relative change in the solution
%           is less than tol
%
%   [X, R, RES] = fitcircle(...) returns the 2 norm of the residual from 
%   the least squares fit
%
%   Example:
%       x = [1 2 5 7 9 3; 7 6 8 7 5 7];
%       % Get linear least squares fit
%       [zl, rl] = fitcircle(x, 'linear')
%       % Get true best fit
%       [z, rad] = fitcircle(x)
%
%   Reference: "Least-squares fitting of circles and ellipses", W. Gander,
%   G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

% This implementation copyright Richard Brown, 2007, but is freely
% available to copy, use, or modify as long as this line is maintained

error(nargchk(1, 5, nargin, 'struct'))

% Default parameters for Gauss Newton minimisation
params.maxits = 1000;
params.maxfuneval = 1000;
params.tol    = 1e-5;
params.maxdiameter = inf;
params.maxxdisp= inf;

% Check x and get user supplied parameters
[x, fNonlinear, params] = parseinputs(x, params, varargin{:});
fNonlinear=true;

% Convenience variables
m  = size(x, 2);
x1 = x(1, :)';
x2 = x(2, :)';
x20 =x2;

[z0 rad0]= simpleCircleGuess(x);


% 2) Nonlinear refinement to miminise geometric error, and compute residual
[z, rad] = fitEllipseFixed(x, z0, rad0,errorRatioYX)
hold off;
plot(x(1,:),x(2,:),'o');
hold all;
plotellipse(z,r,r,0,'k-');
axis equal

% END MAIN FUNCTION BODY

% NESTED FUNCTIONS
    function [z, rad] = fitEllipseFixed(x, z0, r0,errorRatioYX)
      maxd = params.maxdiameter;

      % Get initial phase estimates
      m = size(x, 2);
      maxd = params.maxdiameter;
      phi0 = angle( [1 i] *(x - repmat(z0, 1, m)) )';
      u0 = [phi0; z0; r0];
      lb=[-inf*ones(size(phi0));[-Inf;-Inf];0];
      ub=[+inf*ones(size(phi0));[Inf;Inf];maxd/2];
        
      %options = optimset('MaxIter',params.maxits,'TolFun',params.tol,'MaxFunEvals',params.maxfuneval,'Display','final');
      options = optimset('MaxIter',params.maxits,'TolFun',params.tol,'MaxFunEvals',params.maxfuneval,'Display','final','Jacobian','On');
      u = lsqnonlin(@sys,u0,lb,ub,options);
      z = u(end-2:end-1);
      rad= u(end);
        
        function [f,J] = sys(u)
            %SYS   Nonlinear system to be minimised - the objective
            %function is the distance to each point from the fitted circle
            %contained in u

            phi = u(1:end-3);
            z = u(end-2:end-1);
            r= u(end);

            % Convenience trig variables
            c = cos(phi);
            s = sin(phi);

            % Objective function  - apply error weighting to the axes
            f = zeros(2 * m, 1);
            J = zeros(2 * m, m + 3);
            for ii =1:m
              f(2*ii-1)= x(1, ii) - z(1) - r* c(ii);
              %f(2*ii)  = (x(2, ii) - z(2) - r* s(ii));
              f(2*ii)  = (x(2, ii) - z(2) - r* s(ii))/errorRatioYX^2;

              % Jacobian
              %rows = [2*ii-1,2*ii];
              %J(rows, ii) = - r * [-1 * s(ii);  c(ii)];
              %J(rows, (end-2:end)) = ...
              %    [[-1;0], [0;-1], -1*[c(ii); s(ii)]];
              J((2*ii-1), ii) = - r * [-1 * s(ii)];
              J((2*ii), ii)   = - r * [ c(ii)]/errorRatioYX^2;
              J((2*ii-1), (end-2:end)) = [-1, 0, -1*c(ii)];
              J((2*ii), (end-2:end)) =   [0, -1, -1*s(ii)]/errorRatioYX^2;
            end
            
        end % sys
        
    end % fitcircle_geometric

% END NESTED FUNCTIONS

end % fitcircle

function [x, fNonlinear, params] = parseinputs(x, params, varargin)
% Make sure x is 2xN where N > 3
if size(x, 2) == 2
    x = x';
end

if size(x, 1) ~= 2
    error('fitcircle:InvalidDimension', ...
        'Input matrix must be two dimensional')
end

if size(x, 2) < 3
    error('fitcircle:InsufficientPoints', ...
        'At least 3 points required to compute fit')
end

% determine whether we are measuring geometric error (nonlinear), or
% algebraic error (linear)
fNonlinear = true;
switch length(varargin)
    % No arguments means a nonlinear least squares with defaul parameters
    case 0
        return
       
    otherwise
        if rem(length(varargin), 2) ~= 0
            error('fitcircle:propertyValueNotPair', ...
                'Additional arguments must take the form of Property/Value pairs');
        end

        % Cell array of valid property names
        properties = {'maxits', 'tol','maxdiameter','maxfuneval','maxxdisp'};

        while length(varargin) ~= 0
            property = varargin{1};
            value    = varargin{2};
            
            % If the property has been supplied in a shortened form, lengthen it
            iProperty = find(strncmpi(property, properties, length(property)));
            if isempty(iProperty)
                error('fitcircle:UnkownProperty', 'Unknown Property');
            elseif length(iProperty) > 1
                error('fitcircle:AmbiguousProperty', ...
                    'Supplied shortened property name is ambiguous');
            end
            
            % Expand property to its full name
            property = properties{iProperty};
            
            switch property
                case 'maxits'
                    if value <= 0
                        error('fitcircle:InvalidMaxits', ...
                            'maxits must be an integer greater than 0')
                    end
                    params.maxits = value;
                case 'tol'
                    if value <= 0
                        error('fitcircle:InvalidTol', ...
                            'tol must be a positive real number')
                    end
                    params.tol = value;
                case 'maxdiameter'
                    if ~isnumeric(value) || value <= 0
                        error('fitcircle:InvalidMaxdiameter', ...
                            'maxdiameter must be a positive real number')
                    end
                    params.maxdiameter= value;
                case 'maxxdisp'
                    if ~isnumeric(value) || value <= 0
                        error('fitcircle:InvalidMaxxdisp', ...
                            'maxxdisp must be a positive real number')
                    end
                    params.maxxdisp= value;
                case 'maxfuneval'
                    if ~isnumeric(value) || value <= 0
                        error('fitcircle:InvalidMaxfuneval', ...
                            'maxfuneval must be an integer greater than 0')
                    end
                    params.maxfuneval = value;
                end % switch property
            varargin(1:2) = [];
        end % while

  end % switch length(varargin)

end %parseinputs

%-------------------------------------------------
function [z0 r]= simpleCircleGuess(x);

ux = mean(x(1,:));
uy = mean(x(2,:));
sx = std(x(1,:));
sy = std(x(2,:));

z0=[ux;uy];
r=sqrt(sx^2+sy^2);
end

