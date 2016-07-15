function generateTformFromData(xy,phPathPre,srPixSize,phPixSize,tformFname,imParamPh,imParamFl)

fprintf('Getting tform points\n');
fprintf('Click the localization first, then the corresponding\n');
fprintf('image position.\n');
fprintf('Press q to stop taking point for this image\n');
fprintf('Press Return to take another point pair\n');

phIm = imread(phPathPre);
phIm = imadjust(phIm, stretchlim(phIm,0),[0 1]);%invert the phase contrast image
hold off;
imshow(phIm);
hold all;
plot(xy(:,1),xy(:,2),'rx','MarkerSize',2);

locPos=[];
imPos=[];
doGetMorePts= true;
while doGetMorePts
  %load the data, plot the tformed overlay
  hLoc = impoint(gca);
  hIm  = impoint(gca);
  str = input('''q'' to stop adding point to this image, Return to add more:','s');
  if strcmp(str,'q')
    doGetMorePts=false;
  end

  lc = getPosition(hLoc);
  ic = getPosition(hIm);
  delete(hLoc);
  delete(hIm);
  locPos = [locPos;lc];
  imPos  = [imPos; ic];
  plot(lc(1),lc(2),'gx');
  plot(ic(1),ic(2),'wo');

end

locPos
imPos

locPosAll={};
imPosAll={};
extraTformHackName = '150206_extraTformHackData.mat';
if exist(extraTformHackName,'file')
  load(extraTformHackName);
end
%append the new localizations
locPosAll = {locPosAll{:},locPos};
imPosAll = {imPosAll{:},imPos};

save(extraTformHackName,'locPosAll','imPosAll');
