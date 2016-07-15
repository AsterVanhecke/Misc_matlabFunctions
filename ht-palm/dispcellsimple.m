function h1=dispcellsimple(cellList,im,varargin)

  nArg = numel(varargin)
  ii=1;
  pixSize=1;
  doLabelCell = false;
  plotArg={};
  while ii<=nArg
    if strcmp(varargin{ii},'PixelSize')
      pixSize=varargin{ii+1};
      ii=ii+2;
    elseif strcmp(varargin{ii},'CellLabel')
      doLabelCell=true;
      cellLabel=varargin{ii+1};
      ii=ii+2;
    else
      plotArg = {plotArg{:},varargin{ii}};
      ii=ii+1;
    end
  end

  h1=gcf;
  hold off;
  imshow(im);
  hold all;
  nCell = numel(cellList{1});
  for ii =1:nCell
    try
      mesh=cellList{1}{ii}.mesh;
      mesh= mesh/pixSize;
      plot(mesh(:,1),mesh(:,2),plotArg{:});
      plot(mesh(:,3),mesh(:,4),plotArg{:});
      if doLabelCell==true
        meshAll = [mesh(:,1:2);mesh(:,3:4)];
        cellCtr = [mean(meshAll,1)];
        text(cellCtr(1),cellCtr(2),['\color{blue}',num2str(cellLabel(ii))]);
      end
    catch ME
     fprintf(  '\n******* Plotting cell %d failed for following reason: *******\n',ii);
     fprintf('%s\n',getReport(ME));
     fprintf(  '**************************************************\n');
   end
  end
  hold off;

