function  plotBwRing(srInMesh,movieName,ringStatName,pixSize,bwName)

minPix = 5;

info = imfinfo(movieName);
nFrame = numel(info);
for ii = 1:nFrame
   im0(:,:,ii) = imread(movieName,'Index',ii,'Info',info);
end

%normalise im0
norm = @(x) (x - min(x(:)))./(max(x(:)) - min(x(:)));

im0 = double(im0);
im0 = norm(im0);

%max individually segmented image too
imBW = zeros(size(im0));%initialise memory
for ii = 1:nFrame
   imC = im0(:,:,ii);
   imBW(:,:,ii) = imC>graythresh(imC);%assign
end

for ii = 1:nFrame
  if ii==1 
     writeMode = 'overwrite';
  else
     writeMode = 'append';
  end
  imwrite(im2uint8(imBW(:,:,ii)),bwName,'WriteMode',writeMode);
  %pause
end
