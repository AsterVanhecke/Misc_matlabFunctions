function BWcrop = cropToRoi(BW,imRoi)
% function BWcrop = cropToRoi(BW,imRoi)
%   imRoi = [ xstart, ystart, xwidth, yheight]

%in case imRoi is outside imageBounds
[sizey,sizex] = size(BW);
imRoi = [max(imRoi(1),1), max(imRoi(2),1), min( imRoi(3),sizex), min(imRoi(4),sizey)];

%crop
crp = BW(imRoi(2):imRoi(2)+imRoi(4)-1, imRoi(1):imRoi(1)+imRoi(3)-1);
% remove edge connected objects
crp = imclearborder(crp);
%return image in original size
BW = BW*0;
BW(imRoi(2):imRoi(2)+imRoi(4)-1, imRoi(1):imRoi(1)+imRoi(3)-1) = crp;
BWcrop = BW;



