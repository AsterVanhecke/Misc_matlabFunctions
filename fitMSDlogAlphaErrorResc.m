function [p,resnorm,residual,exitflag,output] = fitMSDlogAlphaErrorResc(t,msd,tE,p0)
% Function to fit MSD to the function form Backlund et al., which takes
% into account super/subdiffusion (parameter alpha) and static and dynamic
% error. 
% INPUT:
% mmsd: 
% tE: exposure time
% OUTPUT: 
% formula: fitted function formula(p,t)
% (from lsqnonlin)
% p: fitted parameters: [D, alpha, s]
%   p(1), 'D': apparent diffusion coefficient
%   p(2), 'alpha': exponent
%   p(3), 's': localization error (standard deviation of position)
% resnorm,residual,exitflag,output: see documentation lsqnonlin

if ~exist('tE','var')
    tE=0.05;
end
% p = [D, alpha, s]
formula=@(p,t) 0.002*p(1)./((p(2)+2).*(p(2)+1).*tE.^2).*((t+tE).^(p(2)+2)+(t-tE).^(p(2)+2)-2*t.^(p(2)+2))-0.004.*p(1).*tE.^p(2)./((p(2)+2).*(p(2)+1))+2.*(p(3)./10).^2;

fun = @(p)(log(msd)-log(formula(p,t))).^2;
fitopts.FunctionTolerance=10^-16;
fitopts.Display='off';
if ~exist('p0','var')
    p0 = [0.4, 0.4, 0.2];
else
    p0(1)=p0(1)/1000;
    p0(3)=p0(3)/10;
end
lb =  [0    , 0  , 0   ];
ub =  [10 , 2  , 10 ];
[p,resnorm,residual,exitflag,output]= lsqnonlin(fun,p0,lb,ub,fitopts);
p(1)=p(1)/1000;
p(3)=p(3)/10;
end