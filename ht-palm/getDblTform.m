function dblTform_InToBase = getDblTform(In,Base,tformType) 

if ~exist('tformType','var')
   tformType = 'similarity';
end

dblTform_InToBase{1} = cp2tform(In,Base,'similarity');
[XBase1,YBase1] =tforminv(dblTform_InToBase{1},Base(:,1),Base(:,2));
dblTform_InToBase{2} = cp2tform(In,[XBase1,YBase1],tformType);

