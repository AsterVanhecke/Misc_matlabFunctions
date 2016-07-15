function [imOut] = preProcessPhIm(im,szObj,wienerBlockSz)

if ~exist('szObj','var')
   szObj= 150;
end
if ~exist('wienerBlockSz','var')
   wienerBlockSz = [5 5];
end


im = im2double(im);
im = imadjust(im, stretchlim(im,0),[1 0]);%invert the phase contrast image

im = bpass(im,'LObject',szObj,'Normalize');

im = wiener2(im,wienerBlockSz);

%im=adapthisteq(im,'NBins',65536,'ClipLimit',0.01,'NumTiles',[16 16],'Distribution','rayleigh','Alpha',0.4);

imOut = imadjust(im, stretchlim(im,0),[1 0]);%invert the phase contrast image

