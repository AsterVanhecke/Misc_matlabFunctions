function demoBW_Estimation1D()
% this is a demo of the 1D bandwidth estimator

% install path to some tools  
pth = './tools' ; rmpath(pth) ; addpath(pth) ;
pth = [pwd, '/drawTools' ] ; rmpath(pth) ; addpath(pth) ;

% let's check if your bandwidth is properly compiled...
checkcompiledBWestimator() ;

% generate some data
% Note: this is the same way the data is generated in another kde also
% available from Mathworks at http://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator
% That KDE tends to underestimate the bandwidth and is slightly slower.
x = [randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55]';
N = size(x,2) ; d = 1 ;

% construct the input for BW estimation
pdf.Mu = x ;
pdf.Cov{N} = {} ;
for i = 1 : N
    pdf.Cov{i} = zeros(d,d) ;
end
pdf.w = ones(1,N)/N ;

% estimate the bandwidth
tic
H = estimateBandwidth( pdf ) ;
toc

% construct the kde
pdf.Mu = x ;
pdf.Cov{N} = {} ;
for i = 1 : N
    pdf.Cov{i} = H ;
end
pdf.w = ones(1,N)/N ;

% display the kde
figure(1) ;
showDecomposedPdf( pdf ) ;
title('Estimated pdf') ;


