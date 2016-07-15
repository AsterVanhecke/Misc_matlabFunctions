function saveHtpalmVisp(srInMesh,folderName)

  if ~exist(folderName,'dir')
    mkdir(folderName);
  end
  nCell = numel(srInMesh)
  for ii = 1:nCell
    ii
    %cell name
    fname =  sprintf('cell_%03d',ii);
    fnameStub = [folderName,filesep(),fname];

    %data
    srData = srInMesh{ii}.localizations;
    xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
    ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
    zcol = findSRField(srData.parInfo,'semantic','position in sample space in Z');
    fcol = findSRField(srData.parInfo,'semantic','frame number');
    photcol =  findSRField(srData.parInfo,'semantic','emission strength');

    x = srData.data(:,xcol);
    y = srData.data(:,ycol);
    z = srData.data(:,zcol);
    f = srData.data(:,fcol);
    phot = srData.data(:,photcol);

    if ~isempty(x)
      saveVisp(fnameStub,x,y,z,phot,f);
    end
  end

function saveVisp(fnameStub,x,y,z,amp,f)
  fnameout= [fnameStub,'.3d']
  fid =fopen(fnameout,'w');
  dataout = [x(:),y(:),z,amp,f];
  fprintf(fid,'%f\t%f\t%f\t%f\t%u\n',dataout');
  fclose(fid);

