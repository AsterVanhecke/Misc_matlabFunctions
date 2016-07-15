function [z, rad,initguess] = fitCircleScaled(x, errorRatioYX, varargin)
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
params.maxits = 10000;
params.maxfuneval = 10000;
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

%rescale the variables
x2 = x2/errorRatioYX;
x(2,:) = x(2,:)/errorRatioYX;
z0=[z0(1); z0(2)/errorRatioYX];

% 2) Nonlinear refinement to miminise geometric error, and compute residual
[z, rad] = fitEllipseFixed(x, z0, rad0,errorRatioYX)
z(2) = z(2)*errorRatioYX;

hold off; 
plot(x1,x20,'o');
hold all;
plotellipse(z,rad,rad,0,'k-');
axis equal
pause

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
        
      options = optimset('MaxIter',params.maxits,'TolFun',params.tol,'MaxFunEvals',params.maxfuneval,'Display','final');
      u = lsqnonlin(@sys,u0,lb,ub,options);
      rad= u(end);
      z = u(end-2:end-1);
        
        function [f] = sys(u)
            %SYS   Nonlinear system to be minimised - the objective
            %function is the distance to each point from the fitted circle
            %contained in u
            %attempt to implement robust fitting

            phi = u(1:end-3);
            r= u(end);
            z = u(end-2:end-1);

            % Objective function
            f = zeros(2 * m, 1);
            for ii =1:m
              f((2*ii-1):(2*ii))  = x(:, ii) - z - r* [cos(phi(ii));  sin(phi(ii))/errorRatioYX];
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

