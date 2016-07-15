function [tc tg wMax fitOut]  = fitFeingoldForLength(lengthTemp,  w, timeLegthMod)

t = timeLegthMod(lengthTemp);
tc0 = 0;
tg0 = max(t);
wMax0 = 1;
a0 = [tc0,tg0,wMax0];
%lb =[0,0,0];
%ub =[max(t),Inf,2];
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';

a = nlinfit(t,w,@feingoldModel2,a0,opts);
tc = a(1);
tg = a(2);
wMax = a(3);

tFit = unique(sort(t));
wFit = feingoldModel2(a,tFit);
fitOut = [tFit(:),wFit(:)];

%----------------------
function wOut = feingoldModel2(a,t)

tc = a(1);
tg = a(2);
wMax = a(3);
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
wOut(isDiv) =  wMax.*sqrt(1- ((tDiv-tc)./(tg-tc)).^2);
%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOut(isPostDiv)=0;
