function [srInMeshOut, isGoodCell] = deleteBadMeshManual(srInMesh,htPalmSummaryFname)

satVal =0.005;
movieData = load(htPalmSummaryFname);

nCell = numel(srInMesh);

isGoodCell = zeros(nCell,1);
for ii=1:nCell
  ii
  cellData = srInMesh{ii};

  %load the PH im
  fovNo = cellData.fovNo; 
  phIm = imread(movieData.fnameList{fovNo}.phPreImName);

  mesh = cellData.meshPix;
  box = cellData.boxPix;

  hold off;
  hF = gcf;
  phImSat = imadjust(phIm,stretchlim(phIm,[satVal, 1-satVal]),[]);

  imshow(phImSat);

  hold on;
  plot(mesh(:,1),mesh(:,2),'r-','LineWidth',1);
  plot(mesh(:,3),mesh(:,4),'r-','LineWidth',1);
  inputOk = false;
  while ~inputOk
    str = input('Good cell? (''g''/''b''):','s');
    if strcmp(str,'g')
      isGoodCell(ii)  = 1;
      inputOk = true;
    elseif strcmp(str,'b')
      isGoodCell(ii)  = 0;
      inputOk = true;
    end
  end
end

srInMeshOut=srInMesh(find(isGoodCell));
