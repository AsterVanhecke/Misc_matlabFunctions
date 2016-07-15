function demoBW_with_preclustering()

% ----------------------------------------------------------- %
% install path to some tools  
pth = './tools' ; rmpath(pth) ; addpath(pth) ;
pth = [pwd, '/drawTools' ] ; rmpath(pth) ; addpath(pth) ;
% let's check if your bandwidth is properly compiled...
checkcompiledBWestimator() ;
% ----------------------------------------------------------- %

demo_type = 1 ; % 1 or 2
switch demo_type
    case 1
        sclx = 1 ; scly = 1 ;
        [x,y] = meshgrid([0:1:30]*sclx,[0:1:30]*scly) ; data = [x(:)';y(:)'] ;
    case 2        
        x = 0:0.001:1 ; y = sin(x*2) ;  data = [x(:)';y(:)'] ;
end

% construct the input for BW estimation
pdf = precluster( data ) ;

% estimate the bandwidth
H = estimateBandwidth( pdf ) ;

% construct the kdes using the calculated bandwidths
kde = constructKDE(pdf.Mu, pdf.w, H, pdf.Cov) ;

figure(1); clf ; 
% plot the kdes as tabulated distributions
subplot(1,2,1) ;
visualizeKDE('kde', kde, 'tabulated', 0) ; 
I_k = visualizeKDE('kde', kde, 'tabulated', 1) ; 
title('Tabulated pdf: Kristan bw') ;
subplot(1,2,2) ; 
imagesc(I_k); title('Tabulated pdf: Kristan bw') ; 
axis equal ; axis tight ; colormap gray ;

% --------------------------------------------------------------------- %
function pdf = precluster( data )

N = round(0.01*size(data,2)) ;
idx = kmeans(data',N) ;

d = size(data,1)  ;
pdf.Mu = zeros(d, N) ;
pdf.Cov{N} = {} ;
pdf.w = zeros(1,N) ;
for i = 1 : N
    dat = data(:, idx==i) ;       
    pdf.Cov{i} = cov(dat') ;
    pdf.Mu(:,i) = mean(dat,2) ;
    pdf.w(i) = size(dat,2) ;
end
pdf.w = pdf.w/sum(pdf.w) ;

