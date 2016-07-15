function visualizeEstimation()

% ----------------------------------------------------------- %
% install path to some tools  
pth = './tools' ; rmpath(pth) ; addpath(pth) ;
pth = [pwd, '/drawTools' ] ; rmpath(pth) ; addpath(pth) ;
% let's check if your bandwidth is properly compiled...
checkcompiledBWestimator() ;
% ----------------------------------------------------------- %

demo_type = 2 ; % 1 or 2

switch demo_type
    case 1
        sclx = 100 ; scly = 200 ;
        [x,y] = meshgrid([1:2:10]*sclx,[1:2:10]*scly) ; data = [x(:)';y(:)'] ;
    case 2        
        x = 0:0.1:1 ; y = sin(x*2);  data = [x(:)';y(:)'] ;
end
w = ones(1,size(data,2)); w = w / sum(w) ;

% construct the input for BW estimation
N = size(data,2) ; d = size(data,1)  ;
pdf.Mu = data ;
pdf.Cov{N} = {} ;
for i = 1 : N
    pdf.Cov{i} = zeros(d,d) ;
end
pdf.w = ones(1,N)/N ;

% estimate the bandwidth
H = estimateBandwidth( pdf ) ;

% construct the kdes using the calculated bandwidths
kde = constructKDE(data, w, H) ;

figure(1); clf ; 
% plot the kdes as tabulated distributions
subplot(1,2,1) ; visualizeKDE('kde', kde, 'tabulated', 0) ; I_k = visualizeKDE('kde', kde, 'tabulated', 1) ; title('Tabulated pdf: Kristan bw') ;
subplot(1,2,2) ; imagesc(I_k); title('Tabulated pdf: Kristan bw') ; axis equal ; axis tight ; colormap gray ;
