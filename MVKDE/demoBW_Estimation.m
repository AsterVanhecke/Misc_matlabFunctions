function demoBW_Estimation()
% this is a demo of the multivariate bandwidth estimator

% install path to some tools  
pth = './tools' ; rmpath(pth) ; addpath(pth) ;

% let's check if your bandwidth is properly compiled...
checkcompiledBWestimator() ;

% generate some data
d = 3 ; N = 100 ;
x = rand(d,N) ;

% construct the input for BW estimation
pdf.Mu = x ;
pdf.Cov{N} = {} ;
for i = 1 : N
    pdf.Cov{i} = zeros(d,d) ;
end
pdf.w = ones(1,N)/N ;

% estimate the bandwidth
H = estimateBandwidth( pdf ) ;

disp('The estimated bandwidth:')
H


