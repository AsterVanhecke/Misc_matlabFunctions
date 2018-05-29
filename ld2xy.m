function [X,Y]=ld2xy(L,D,mesh)
% Convert cell coordinates L, D (Length, Diameter) to X,Y coordinates
cLine=([mean([mesh(:,1), mesh(:,3)],2) , mean([mesh(:,2), mesh(:,4)],2) ]);
lStep=1; %let's keep this for now
L=L/lStep; % what is the right value?
segCent=cLine(1:end-1,:)+diff(cLine)./2;
segWvect=[mesh(1:end-1,1:2)+diff(mesh(:,1:2))./2 - segCent; 0, 0];
segWvect=segWvect./(segWvect(:,1).^2+segWvect(:,2).^2); %normalize vector
segLengthVec =[diff(cLine); 0, 0];

% precise l
lRem = rem(L,1); % (L , one). l and 1 are too similar l1l1

lIdx=floor(L);
lIdx(lIdx==0)=1;
xLine=cLine(lIdx,1) + (lRem(:)-0.5).*segLengthVec(lIdx,1);
yLine=cLine(lIdx,2) + (lRem(:)-0.5).*segLengthVec(lIdx,2);

X=xLine+D(:).*segWvect(lIdx,1);
Y=yLine+D(:).*segWvect(lIdx,2);
X=reshape(X,size(L));
Y=reshape(Y,size(L));

