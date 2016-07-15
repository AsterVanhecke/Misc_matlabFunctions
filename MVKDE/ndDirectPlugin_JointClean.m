%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2010 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2010
%%
function [H, h_amise, alpha_scale] = ndDirectPlugin_JointClean( Mu, Cov, w, Cov_smp, N_eff )                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: approximation works best when the pdf is spherised ! F = eye() !
d = size(Mu,1) ;
% [new_mu, Cov_smp, w_out] = momentMatchPdf(pdf.Mu, pdf.ps.Cov, pdf.w) ;

% test
% for i = 1 : length(Cov)
%     Cov{i} = diag(diag(Cov{i})) ;
% end
%end test
 
G = (Cov_smp *(4/((d+2)*N_eff))^(2/(d+4))) ; %*0.8^2 ;
% alpha_scale = 1 / det(Cov_smp)^(1/d)  ; 
alpha_scale = 1 ;
F = Cov_smp * alpha_scale; % for numerical stability. it could have been: F = Cov_smp ;
% could also constrain to say that F = identity!
applyapproximation = 1;

if applyapproximation == 1
    % uses approximation !
    % test
    Ccov = zeros(size(Mu,1),size(Mu,2)) ;
    for i = 1 : length(Cov)
        Ccov(:,i) = diag(Cov{i})' ;
    end
    g = (4/((d+2)*N_eff))^(2/(d+4)) ;
    Rf2 = mex_getIntSquaredHessian(double(Mu), double(w), double(Ccov), double(g)) ;
else
    Rf2 = getIntSquaredHessian( Mu, w, Cov, F, G ) ;
end
% if abs(Rf2-R_x) > 1e-3
%     df = 34 ;
% end
% [Rf2, R_x]
%end test

h_amise = (N_eff^(-1) *det(F)^(-1/2) /( sqrt(4*pi)^d * Rf2 * d ))^(1/(d+4)) ;
H = getBWfromStruct( F, h_amise)*alpha_scale ; 
if nargout == 1
    alpha_scale = [] ;
    h_amise = [] ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 





