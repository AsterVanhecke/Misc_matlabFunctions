function [t,w,d,cellData]=pltMTS_diameter2(srInMesh,varargin)

manClass=[];
doBadCell=false;
doManualClassify = false;
doTMax=false;
doFilterList=false;
tMax=Inf;
nArg=numel(varargin);
ii = 1;
while ii <= nArg
  if strcmp(varargin{ii},'ManualClassification')
    doManualClassify=true;
    manClass= varargin{ii+1};
    ii=ii+2;
  elseif strcmp(varargin{ii},'TMax')
    doTMax = true;
    tMax= varargin{ii+1};
    ii=ii+2;
  elseif strcmp(varargin{ii},'FilterList')
    doFilterList= true;
    filterList= varargin{ii+1};
    ii=ii+2;
  else
    ii=ii+1;
  end
end

nCell = numel(srInMesh);
kk = 1;
for ii =1:nCell

  tExpt(kk) = srInMesh{ii}.tExpt_min;
  medD = srInMesh{ii}.diameter.medD;
  conSiteD = srInMesh{ii}.diameter.conSiteD;
  plotMinD(kk) = conSiteD(1,2);
  plotMaxD(kk) = medD;
  W(kk) = plotMinD(kk)/plotMaxD(kk);
  fov(kk) =srInMesh{ii}.fovNo;
  kk=kk+1;
end

cellOk=ones(nCell,1);

if tMax 
  tOk = tExpt<tMax;
  cellOk = cellOk & tOk(:);
end
if doFilterList
  cellOk = cellOk & filterList;
end
if doManualClassify
  %first exclude the bad cell
  isBad=logical(zeros(size(cellOk)));
  isBad(manClass.badCell) = 1;
  cellOk = cellOk & ~isBad;
  %annotate the septated cells as 'bad cells' 
  %and plot separately
  isDiv = logical(zeros(size(cellOk)));
  isDiv(manClass.postDiv)=1;
  %HACK - list the cells including diameter
  cellData = [(1:nCell)',tExpt(:),plotMinD(:),W(:)];
  cellData(isDiv,3:4)=0;
  cellData(~cellOk,:)=[];

  %END HACK

  %separate classification keeping the manual cells
  cellOk = cellOk & ~isDiv;

  %annotate the manual cells
  tDiv = tExpt(isDiv);
  tDiv = tDiv(:);
  dDiv = zeros(size(tDiv));
  dDiv2 = [plotMinD(isDiv)', plotMaxD(isDiv)', W(isDiv)'];
  %add the non-segmented post-div cells
  pdExt = manClass.postDivExtra;
  nFov = size(pdExt,1);
  for ii = 1:nFov
    nExtra= pdExt(ii,3);
    if nExtra>=0
      tExt = pdExt(ii,2)*ones(nExtra,1);
      dExt = zeros(nExtra,1);
      dExt2 = zeros(nExtra,3);
      tDiv = [tDiv;tExt];
      dDiv = [dDiv;dExt];
      dDiv2 = [dDiv2;dExt2];
    end
  end

end

figure;
subplot(2,1,1);
hold on;
plot(tExpt(cellOk),plotMinD(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv,'ko');
end
ylim([0 600]);

xlabel('T (min)');
ylabel('Constriction site diameter (nm)');
subplot(2,1,2);
hold on;
plot(tExpt(cellOk),W(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv,'ko');
end
ylim([0 1.2]);
ylabel('Waist ratio r/R');
xlabel('T (min)');

figure;
subplot(2,1,1);
hold on;
plot(tExpt(cellOk),plotMinD(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv2(:,1),'ro');
end
ylim([0 600]);

xlabel('T (min)');
ylabel('Min. diameter (nm)');
subplot(2,1,2);
hold on;
plot(tExpt(cellOk),W(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv2(:,3),'ro');
end
ylim([0 1.2]);
ylabel('Waist ratio r/R');
xlabel('T (min)');


%%try a contour plot
%t= [tExpt(cellOk)';tDiv(:)];
%w= [W(cellOk)';dDiv(:)];
%tCtr = 0:20:max(t);
%wCtr = 0:0.05:1.2;
%ctrs={tCtr,wCtr}
%n=hist3([t,w],'Ctrs',ctrs);
%NMAX = 20;
%n(n>NMAX)=NMAX;
%figure;
%nCountour=16;
%h = contourf(tCtr,wCtr,n');
%colormap('gray');
%cmap = colormap;
%cmap = flipud(cmap);
%colormap(cmap);
%ylim([0 1.2]);
%ylabel('Waist ratio r/R');
%xlabel('T (min)');
%
%figure;
%h = pcolor(tCtr,wCtr,n');
%set(h,'edgecolor','none');
%colormap('gray');
%cmap = colormap;
%cmap = flipud(cmap);
%colormap(cmap);
%ylim([0 1.2]);
%ylabel('Waist ratio r/R');
%xlabel('T (min)');
%

%try an area plot
t= [tExpt(cellOk)';tDiv(:)];
w= [W(cellOk)';dDiv(:)];
d = [plotMinD(cellOk)';dDiv(:)];
a = pi*d.^2/4;

figure;
plot(t,a,'ko');
xlabel('T (min)');
ylabel('Area (nm^2)');

