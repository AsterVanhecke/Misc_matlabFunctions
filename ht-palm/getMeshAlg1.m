function cellList = getMeshAlg1(img,regions0,p)

  %initialize cellList
  cellList = cell(1);
  cellList{1} = cell(1);

  %initialize strel
  [se maskdx maskdy] = initP();

  imsizes=[size(img,1) size(img,2)];

  if p.invertimage
      img = max(img(:)) - img;
  end

  %adapted from processFrameOld
  imge = img2imge(img,p.erodeNum);
  imge16 = img2imge16(img,p.erodeNum);
  imgemax = max(max(imge));
  img = img2imge(img,0);




  regions0 = bwlabel(regions0,4);

  stat0 = regionprops(regions0,'area','boundingbox');
  reg=1;
  cell0 = 0;
  frame = 1;
  regmax = length(stat0);
  regold = 0;
  regmindisp=true;

  while reg<=regmax 
    if reg>regold, repcount=0; else repcount=repcount+1; end
    if repcount>20, reg=reg+1; continue; end
    regold = reg;
    if regmindisp, display(['processing region ' num2str(reg)]); end
    if reg>length(stat0), stat0=regionprops(regions0,'area','boundingbox'); end
    if reg>length(stat0), break; end
    statC = stat0(reg);
    
    % otherwise compute the properties for proper splitting
    roiBox(1:2) = max(statC.BoundingBox(1:2)-p.roiBorder,1); % coordinates of the ROI box
    roiBox(3:4) = min(statC.BoundingBox(1:2)+statC.BoundingBox(3:4)+p.roiBorder,[size(img,2) size(img,1)])-roiBox(1:2);
    roiRegs = imcrop(regions0,roiBox); % ROI with all regions labeled
    roiMask = roiRegs==reg; % ROI with region #reg labeled
    % for k=1:p.dilMaskNum, roiMask = imdilate(roiMask,se); end
    roiImg = imcrop(imge,roiBox);
    roiImgI = imcrop(img,roiBox);
    
    perim = bwperim(imdilate(roiMask,se));
    pmap = 1 - perim;
    f1=true;
    while f1
      pmap1 = 1 - perim + imerode(pmap,se);
      f1 = max(max(pmap1-pmap))>0;
      pmap = pmap1;
    end;
    regDmap = (3*(pmap+1) + 5*(1-roiImg/imgemax)).*roiMask;
    maxdmap = max(max(regDmap));


    % Get and smooth the boundary of the region
    [ii,jj]=find(bwperim(roiMask),1,'first');
    cCell=bwtraceboundary(roiMask,[ii,jj],'n',8,inf,'counterclockwise');
    if isfield(p,'interpoutline') && p.interpoutline
      ctrlist = interpoutline(cCell,roiImg,p);
    else
      ctrlist = {cCell};
    end

    if isempty(ctrlist), display(['fitting region ' num2str(reg) ' - failed, no outline creating during smoothing']); end
    for cind = 1:length(ctrlist)
      cCell = ctrlist{cind};
      if p.fsmooth>0 && p.fsmooth<Inf
        fpp = frdescp(cCell);
        cCell = ifdescp(fpp,p.fsmooth);
        cCell=[cCell;cCell(1,:)]; % make the last point the same as the first one
      end
      cCell = rot90(cCell,2);
      cCell(:,1) = cCell(:,1)+roiBox(1)-1;
      cCell(:,2) = cCell(:,2)+roiBox(2)-1;
      
      % Save the data
      if p.getmesh, mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth); end

      if (~p.getmesh || length(mesh)>4) && ...
         min(cCell(:,1))>= p.imRoi(1) &&...
         min(cCell(:,2))>= p.imRoi(2) &&...
         max(cCell(:,1))<=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)) && ...
         max(cCell(:,2))<=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)) 
        cellStruct.algorithm = 1;
        if p.getmesh, cellStruct.mesh = mesh; end
        cellStruct.polarity = 0;
        cellStruct.box = roiBox;
        cellStruct.ancestors = [];
        cellStruct.descendants = [];
        cellStruct.divisions = [];
        cellStruct.stage = [];
        cellStruct.contour = cCell;
        cellStruct.polarity = 0; % algorithm 1 cannot divide cells unless this is interpoutline

        cellStruct.birthframe = frame;
        cellStruct = getextradata(cellStruct);%add length area etc info

        cell0 = cell0+1;
        cellList{frame}{cell0} = cellStruct;
        display(['fitting region ' num2str(reg) ' - passed and saved as cell0 ' num2str(cell0)])
      else
        % if the cell0 is not on the image border OR it did not pass the
        % quality test - split it
        reason = 'unknown';
        if p.getmesh && length(mesh)<=4, reason = 'bad region shape quality'; end
        if min(cCell(:,1))<=p.imRoi(1), reason = 'cell0 on x=0 boundary'; end
        if min(cCell(:,2))<=p.imRoi(2), reason = 'cell0 on y=0 boundary'; end
        if max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)), reason = 'cell0 on x=max boundary'; end
        if max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)), reason = 'cell0 on y=max boundary'; end
        display(['fitting region ' num2str(reg) ' - quality check failed - ' reason])
      end
    end
    reg=reg+1;
  end


end
%------------------------------------------------------
function [se maskdx maskdy] = initP()
  % initModel
  
  se = strel('arb',[0 1 0;1 1 1;0 1 0]); % erosion mask, can be 4 or 8 neighbors
  maskdx = fliplr(fspecial('sobel')'); %[-1 0 1; -2 0 2; -1 0 1]/2; % masks for computing x & y derivatives
  maskdy = fspecial('sobel');%[1 2 1; 0 0 0; -1 -2 -1]/2;
end

%-----------------------------------------
%-----------------------------------------

function im2 = img2imge16(im,nm)
% erodes image "im" by "nm" pixels
se = strel('arb',[0 1 0;1 1 1;0 1 0]); % erosion mask, can be 4 or 8 neighbors
im2 = im;
for i=1:nm, im2 = imerode(im2,se);end
end

%-----------------------------------------
%-----------------------------------------


function im2 = img2imge(im,nm)
% erodes image "im" by "nm" pixels and normalizes it, outputs double
se = strel('arb',[0 1 0;1 1 1;0 1 0]); % erosion mask, can be 4 or 8 neighbors
im2 = im;
for i=1:nm, im2 = imerode(im2,se);end
im2 = double(im2);
mn=mmin(im2);
mx=mmax(im2);
im2=1-(im2-mn)/double(mx-mn);
end


function [out]= mmin(array)
    %
    % returns the min value of a N dimensional array
    %
    out=min(array);
    n=ndims(array);
    if n>1
        for i=2:n
            out=min(out);
        end
    end
end

function [out]= mmax(array)
    %
    % returns the max value of a N dimensional array
    %
    out=max(array);
    n=ndims(array);
    if n>1
        for i=2:n
            out=max(out);
        end
    end
end

%% Cell statistics and additional data

function str=getextradata(str)
    % calculates geometrical properties of detected meshes in structure
    % "str", corresponding to a single cell in "cellList", recording to the
    % same structure. The properties are: "steplength", "steparea",
    % "stepvolume" (length, area, and volume of each segment), "length",
    % "area", "volume" (foe the whole cell), "lengthvector" - coordinates
    % along cell centerline of the centers of each segment.
    if isempty(str), return; end
    if isfield(str,'mesh') && length(str.mesh)>1
        mesh = str.mesh;
        lng = size(mesh,1)-1;
        if ~isfield(str,'polarity'), str.polarity=0; end

        %length
        str.steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
            mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
        str.length = sum(str.steplength);

        str.lengthvector = cumsum(str.steplength)-str.steplength/2;

        % area
        str.steparea = [];
        for i=1:lng, str.steparea=[str.steparea;polyarea([mesh(i:i+1,1);mesh(i+1:-1:i,3)],[mesh(i:i+1,2);mesh(i+1:-1:i,4)])]; end
        str.area = sum(str.steparea);

        % volume
        d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
        str.stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*str.steplength*pi/4;
        str.volume = sum(str.stepvolume);
    elseif isfield(str,'contour') && length(str.contour)>1
        contour = str.contour;
        lng = size(contour,1);
        str.length = sqrt(max(max( (repmat(contour(:,1),1,lng)-repmat(contour(:,1)',lng,1)).^2 + ...
                                   (repmat(contour(:,2),1,lng)-repmat(contour(:,2)',lng,1)).^2)));
        str.area = polyarea(contour(:,1),contour(:,2));
    end
end

function d=edist(x1,y1,x2,y2)
    % complementary for "getextradata", computes the length between 2 points
    d=sqrt((x2-x1).^2+(y2-y1).^2);
end

