%-------------------------------------------------------------------------%
%                main function: fitAreaSeamusModBound                         %
%-------------------------------------------------------------------------%
% This function allows to compute the fit of the waist diamtert w versus
% time considering a model where
% - the area inserted at septum is a 2D circle
% - the area rate insertion is constant
%
% The input of the function are:
% - w^2 (the square of the diameter at septum divided by the max diameter)
% - t (time)
% The output of the function are the fitting parameters

function [tc tg wMax fitOut]  = fitAreaSeamusModBound(t, w2)

tc0 = 0;
tg0 = max(t);
wMax0 = 1;
a0 = [tc0, tg0, wMax0];
lb =[0,0,0];
ub =[max(t),Inf,1];
a = lsqcurvefit(@seamusAreaModel, a0, t, w2, lb, ub);
tc = a(1);
tg = a(2);
wMax = a(3);

tFit = unique(sort(t));
wFit = seamusAreaModel(a, tFit);
fitOut = [tFit(:), wFit(:)];
end

%----------------------
function wOutSquared = seamusAreaModel(a, t)

tc = a(1);
tg = a(2);
wMax = a(3);
%fic

%3 conditions
%1. (t<tc), w=wMax
%2. (tc<t<tg), w=linear constriction
%3. (t>tg) w = 0;

wOutSquared = 0.*t;

%1. (t<tc), w=wMax
isPreDiv = t<tc;
wOutSquared(isPreDiv) = wMax^2;
%wOut(isPreDiv)=1;

%2. (tc<t<tg),w=linear constriction
%passes through points
%(tc,wMax), (tg,0);
isDiv = (t >= tc & t < tg);
tDiv = t(isDiv);
wOutSquared(isDiv) = wMax^2.*( 1 - 1/(tg-tc).*(tDiv-tc) );


%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOutSquared(isPostDiv)=0;
end

