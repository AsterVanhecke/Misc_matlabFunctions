function [cellListOut]= associatePhSr(cellList, srData);
% put SR points in PH meshes (use spotfinder routines)


nFrame = numel(cellList);
%for each frame
for ii = 1:nFrame
  nCell = numel(cellList{ii});
  %for each cell
  for jj = 1:nCell
    if ~isempty(cellList{ii}{jj})
      box= cellList{ii}{jj}.box;
%%<DEBUG
%      hold off;
%      mesh = cellList{ii}{jj}.mesh;
%      %mesh = mesh/pixSize;
%      plot(mesh(:,1),mesh(:,2),'g-')%left side of cell
%      hold all;
%      plot(mesh(:,3),mesh(:,4),'g-');%right side of cell
%      boxV = box2Vertex(box);
%      plot([boxV(:,1);boxV(1,1)],[boxV(:,2);boxV(1,2)],'r-');
%%DEBUG>
      %get SR points within the box
      srDataBox = findInBox(srData,box);
      %parse within box points and add to cellList structure
      cellList{ii}{jj} = addSr2Cell(cellList{ii}{jj},srDataBox);

    end
  end
end

cellListOut = cellList;

%-------------------------------------------------------
function srDataBox = findInBox(srData,box)
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');

rowInBox = find(srData.data(:,xcol)>box(1) & srData.data(:,xcol)<(box(1)+box(3)-1) & srData.data(:,ycol)>box(2) & srData.data(:,ycol)<(box(2)+box(4)-1)  );

srDataBox.parInfo = srData.parInfo;
srDataBox.data    = srData.data(rowInBox,:);

%%<DEBUG
%hold off;
%plot(srData.data(:,xcol),srData.data(:,ycol),'kx');
%hold on;
%plot(srDataBox.data(:,xcol),srDataBox.data(:,ycol),'gx');
%boxV = box2Vertex(box);
%plot([boxV(:,1);boxV(1,1)],[boxV(:,2);boxV(1,2)],'r-');
%pause;
%%DEBUG>


%-------------------------------------------------------
function cellOut= addSr2Cell(cellIn,srData)

m=mTrack;
xcol = findSRField(srData.parInfo,'semantic','position in sample space in X');
ycol = findSRField(srData.parInfo,'semantic','position in sample space in Y');
   
cellIn.localizations = srData;

%convert cell outline to cell coordinate system
cellIn.meshAligned = cellIn.mesh;
[l1 d1] = projectToMeshMod(cellIn.mesh(:,1),cellIn.mesh(:,2),cellIn.mesh,cellIn.steplength);%get signed cell coordinates
[l2 d2] = projectToMeshMod(cellIn.mesh(:,3),cellIn.mesh(:,4),cellIn.mesh,cellIn.steplength);%get signed cell coordinates
cellIn.meshAligned(:,1) =l1;
cellIn.meshAligned(:,2) =d1;
cellIn.meshAligned(:,3) =l2;
cellIn.meshAligned(:,4) =d2;

%calculate position in cell coordinates
X = [srData.data(:,xcol),srData.data(:,ycol)];
[cellIn.l cellIn.d] = projectToMeshMod(X(:,1),X(:,2),cellIn.mesh,cellIn.steplength);%get signed cell coordinates

%calculate segment number
%have to change to box coordinate system
box =  cellIn.box;
XBox = ...
  [(X(:,1) - box(1) +1), (X(:,2) - box(2) +1) ];%box coordinate system is euclidean from top left image corner =(1,1)
   
cellIn.position = m.assignSegment(XBox(:,1),XBox(:,2),cellIn.mesh,box);

%%<DEBUG
%hold off
%mesh = cellIn.meshAligned;
%plot(mesh(:,1),mesh(:,2),'g-')%left side of cell
%hold all;
%plot(mesh(:,3),mesh(:,4),'g-');%right side of cell
%IN = cellIn.position>0;
%%IN = cellIn.position>-1;
%%plot(cellIn.l(IN),cellIn.d(IN),'kx');
%plot(cellIn.l,cellIn.d,'kx');
%axis equal;
%pause
%%DEBUG>
cellOut=cellIn;
