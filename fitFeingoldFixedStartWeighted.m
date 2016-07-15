function [tg wMax fitOut]  = fitFeingoldFixedStartWeighted(t, w, tc, weights)

tg0 = max(t);
wMax0 = 1;
a0 = [tg0,wMax0];

opts = statset('nlinfit');
%opts.RobustWgtFun = 'bisquare';

a = nlinfit(t,w,@(x,xdata) feingoldModel2(x,xdata,tc),a0,opts, 'Weights', weights);
tg = a(1);
wMax = a(2);

tFit = unique(sort(t));
wFit = feingoldModel2(a,tFit,tc);
fitOut = [tFit(:),wFit(:)];

%----------------------
function wOut = feingoldModel2(a,t,tc)

tg = a(1);
wMax = a(2);
%fic

%3 conditions
%1. (t<tc), w=wMax
%2. (tc<t<tg), w=feingold model
%3. (t>tg) w = 0;

wOut = 0.*t;

%1. (t<tc), w=wMax
isPreDiv = t<tc;
wOut(isPreDiv)=wMax;

%2. (tc<t<tg), w=feingold model
isDiv = (t>=tc& t<tg);
tDiv = t(isDiv);
wOut(isDiv) = wMax.*sqrt(1- ((tDiv-tc)./(tg-tc)).^2);

%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOut(isPostDiv)=0;
