function [formula, p,resnorm,residual,exitflag,output] = fitMSDlogAlphaError(t,msd,tE,p0)
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

% To prevent t(1)<tE by a tiny amount, leading to imaginary solutions
% because of the term "(t-tE).^(p(2)+2)"
tE=tE-10^-16;

formula=@(p,t) 2*p(1)./((p(2)+2).*(p(2)+1).*tE.^2).*((t+tE).^(p(2)+2)+(t-tE).^(p(2)+2)-2*t.^(p(2)+2))-4.*p(1).*tE.^p(2)./((p(2)+2).*(p(2)+1))+2.*p(3).^2;

fun = @(p)(log(msd)-log(formula(p,t))).^2;
fitopts.FunctionTolerance=10^-16;
fitopts.Display='off';
if ~exist('p0','var')
    p0 = [0.0004, 0.4, 0.02];
end
lb =  [0    , 0  , 0   ];
ub =  [0.001 , 2  , 0.1 ];
[p,resnorm,residual,exitflag,output]= lsqnonlin(fun,p0,lb,ub,fitopts);
end