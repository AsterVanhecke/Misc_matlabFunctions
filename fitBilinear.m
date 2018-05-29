function [A, B, C, tc, fitOut, resnorm, residual] = fitBilinear(t,w)
%% Initial and bounding parameters
A0=0;
B0=nanmedian(w);
C0=0;
tc0=nanmedian(t);
a0=[A0, B0, C0, tc0]; %initial values for fit
% lb=[-Inf, -Inf, -Inf]; %lower bounds for fit
% ub=[Inf, Inf, Inf, Inf]; %upper bounds for fit
%% Fitting
[a, resnorm, residual]= lsqcurvefit(@bilinearModel, a0, t, w);
A=a(1);
B=a(2);
C=a(3);
tc=a(4);

tFit = unique(sort(t));
wFit = bilinearModel(a,tFit);
fitOut = [tFit(:),wFit(:)];
end

function wOut = bilinearModel(a, tIn)
A=a(1);
B=a(2);
C=a(3);
tc=a(4);
wOut=nan(size(tIn));
wOut(tIn<tc)=A.*tIn(tIn<tc)+B;
wOut(tIn>=tc)=C.*tIn(tIn>=tc)+B+(A-C)*tc;
end