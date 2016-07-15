function [tc, tg, wMax, alpha, fitOut, resnorm, residual, BIC]  = fitXiao(varargin)
% function [tc, tg, wMax, alpha, fitOut, resnorm, residual, BIC]  = fitXiao(t, w)
% This function allows to compute the fit of the waist diameter, w, versus
% time, t, to the Jiao model Xiao
% The input of the function are:
% - w (diameter at septum divided by the max diameter)
% - t (time)
% The output of the function are the fitting parameters, tc, tg, wMax and
% alpha, as weel as the points making up the fitted curve. Optional output:
% resnorm and residual of the lsqcurvefit and the Bayesian information
% criterion.
narginchk(2,3);
t=varargin{1,1};
w=varargin{1,2};
if size(varargin,2) == 2
    lAlpha=0.1;
    uAlpha=10;
    nparam=4;
    alpha0 = 2;
elseif size(varargin,2) == 3
    lAlpha=varargin{1,3};
    uAlpha=lAlpha;
    nparam=3;
    alpha0 = lAlpha;
else
    error('Wrong number of input parameters, correct format is fitXiao(timevec, widthvec) or fitXiao(timevec, widthvec, alpha)') 
end
tc0 = 0;
tg0 = max(t);
wMax0 = 1;
a0 = [tc0,tg0,wMax0, alpha0];
lb =[0,0,0.95,lAlpha];
ub =[max(t),Inf,1.2,uAlpha];
[a, resnorm, residual] = lsqcurvefit(@XiaoModel, a0, t, w, lb, ub);
tc = a(1);
tg = a(2);
wMax = a(3);
alpha = a(4);

% Compute the Bayesian information criterion
numDataPoints = length(residual);
BIC = numDataPoints*log(resnorm)+nparam*log(numDataPoints);

tFit = unique(sort(t));
wFit = XiaoModel(a,tFit);
fitOut = [tFit(:),wFit(:)];
end

%----------------------
function wOut = XiaoModel(a, t)

tc = a(1);
tg = a(2);
wMax = a(3);
alpha = a(4);
%fic

%3 conditions
%1. (t<tc), w=wMax
%2. (tc<t<tg), w=linear constriction
%3. (t>tg) w = 0;
%4. tc<tg, so punish when tg<t<tc

wOut = 0.*t;

%1. (t<tc), w=wMax
isPreDiv = t<tc;
wOut(isPreDiv) = wMax;
%wOut(isPreDiv)=1;

%2. (tc<t<tg),w=linear constriction
%passes through points
%(tc,wMax), (tg,0);
isDiv = (t>=tc & t<tg);
tDiv = t(isDiv);
wOut(isDiv) = wMax.*( 1 - ((tDiv-tc)./(tg-tc)).^alpha ).^(1./alpha);


%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOut(isPostDiv)=0;

%4. tc<tg, so punish when tg<t<tc
isYouDoneFdUp = (t>=tg & t<tc);
wOut(isYouDoneFdUp)=-9001;
end