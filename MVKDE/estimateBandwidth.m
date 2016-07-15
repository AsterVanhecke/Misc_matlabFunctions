% Author: Matej Kristan, Faculty of Computer and Information Science, University of Ljubljana (2013).
% Email: matej.kristan@fri.uni-lj.si
% 
% This is a faster (approximate) bandwidth estimation that follows [1]. The approximation
% does not appear to affect much the accuracy, but it executes faster 
% than its exact counterpart in Matlab implementation. 
% 
% [1] Kristan Matej, Leonardis Ales and Skocaj Danijel, "Multivariate Online Kernel Density Estimation with Gaussian Kernels", 
% Pattern Recognition, 2011.
%
% If you use this code, please cite the paper [1].
%
function H = estimateBandwidth( pdf0, N_eff )
%
% pdf0 ... the input distribution -- note that the bandwidth is usually
% estimated directly from the samples. But sometimes it is beneficial to
% first precluster data, approximate each cluster by a Gaussian and
% estimate the bandwidth from that. The standard bandwidth estimators will
% not allow that, but the bandwidth from Kristan et al., [1] does. If you
% want to estimate the bandwidth directly from the data, just set the
% means in the distribution to your data points and the covariances to
% zero.
% pdf0.Mu  ... a dXN array, with N d-dimensional mean values of the
%              components
% pdf0.Cov ... N-component cell array of covariances, with each cell being dXd array
% pdf0.w   ... 1XN array of component weights
%
% N_eff ... is the effective sample size. If this is left empty, the number
% of components in the pdf will be taken as N_eff. However, if you have
% preclustered the data, then you should set the N_eff to the number of
% original samples in clustering.

if nargin < 2
    N_eff = length(pdf0.w) ;
end

% first we'll spherize the distribution
[Mu, C, ~] = momentMatchPdf(pdf0.Mu, pdf0.Cov, pdf0.w) ;
[U, S, V] = svd(C) ; 
T = (diag(1./sqrt(diag(S))))*U' ;

pdf0 = applyForScaleTransformToPdf( pdf0, Mu, T ) ;
C = T*C*T' ;

% then we'll calculate the optimal bandwidth by Kristan's estimator
[H, ~, ~] = ndDirectPlugin_JointClean( pdf0.Mu, pdf0.Cov, pdf0.w, C, N_eff ) ;

pdf.Mu = Mu ;
pdf.Cov = {H} ;
pdf.w = 1 ;
pdf = applyInvScaleTransformToPdf( pdf, Mu, T ) ; 
H = pdf.Cov{1} ;
