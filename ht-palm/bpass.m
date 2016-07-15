function res= bpass(image_array,varargin)
% 
% NAME:
%               bpass
% PURPOSE:
%               Implements a real-space bandpass filter that suppresses 
%               pixel noise and long-wavelength image variations while 
%               retaining information of a characteristic size.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               res = bpass( image_array, lnoise, lobject )
% INPUTS:
%               image:  The two-dimensional array to be filtered.
%               ..., 'LNoise', lnoise: Characteristic lengthscale of noise in pixels.
%                       Additive noise averaged over this length should
%                       vanish. May assume any positive floating value.
%                       May be set to 0 or false, in which case only the
%                       highpass "background subtraction" operation is 
%                       performed.
%               ..., 'LObject', lobject: (optional) Integer length in pixels somewhat 
%                       larger than a typical object. Can also be set to 
%                       0 or false, in which case only the lowpass 
%                       "blurring" operation defined by lnoise is done,
%                       without the background subtraction defined by
%                       lobject.  Defaults to false.
%                       performed at all. 
%               ..., 'Normalize': normalise image between 0-1, return type double
%
% OUTPUTS:
%               res:    filtered image.
% PROCEDURE:
%               simple convolution yields spatial bandpass filtering.
% NOTES:
% Performs a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 
% original - lowpassed version of original => highpassed version of the
% original
% 
% Performing a lowpass and a highpass results in a bandpassed image.
% 
% Converts input to double.  Be advised that commands like 'image' display 
% double precision arrays differently from UINT8 arrays.
% 
% Based on Crocker and Grier's bpass.m
% SH - added padding so you dont chop off the image edges
% SH -renormalise to 0-1 at end

lobject = false;
lnoise = false;
doNorm = false;
ii = 1;
while ii <= numel(varargin)
   if strcmp(varargin{ii},'LNoise')
      lnoise= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'LObject')
      lobject = varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'Normalize')
      doNorm = true;
      ii = ii+1;
   else
      ii = ii+1;
   end
end

imType = class(image_array);
imIsInt = isinteger(image_array);
lzero = max(lobject,ceil(5*lnoise));
%pad the image - this means that we dont lose the corners when we do the convolution
image_array = padarray(image_array,[lzero lzero],'symmetric');

image_array = double(image_array);

normalize = @(x) x/sum(x);
if lnoise == 0
  gaussian_kernel = 1;
else      
  gaussian_kernel = normalize(...
    exp(-((-ceil(5*lnoise):ceil(5*lnoise))/(2*lnoise)).^2));
end

if lobject  
  boxcar_kernel = normalize(...
      ones(1,length(-round(lobject):round(lobject))));
end
  
% JWM: Do a 2D convolution with the kernels in two steps each.  It is
% possible to do the convolution in only one step per kernel with 
%
  % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
  % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
% 
% but for some reason, this is slow.  The whole operation could be reduced
% to a single step using the associative and distributive properties of
% convolution:
%
  % filtered = conv2(image_array,...
  %   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
  %   'same');
%
% But this is also comparatively slow (though inexplicably faster than the
% above).  It turns out that convolving with a column vector is faster than
% convolving with a row vector, so instead of transposing the kernel, the
% image is transposed twice.
 
if lnoise
   gconv = conv2(image_array',gaussian_kernel','same');
   gconv = conv2(gconv',gaussian_kernel','same');
else
   gconv = image_array;
end

if lobject
  bconv = conv2(image_array',boxcar_kernel','same');
  bconv = conv2(bconv',boxcar_kernel','same');

  filtered = gconv - bconv;
else
  filtered = gconv;
end

% Delete the padded edges now we're finished

filtered((end - lzero + 1):end,:) = [];
filtered(:,(end - lzero + 1):end) = [];
filtered(1:(round(lzero)),:) = [];
filtered(:,1:(round(lzero))) = [];

if doNorm == true
   %normalise to between 0-1
   res= (filtered-min(filtered(:)))/ (max(filtered(:)) - min(filtered(:)));
elseif strcmp(imType,'double')
   res = filtered;
elseif imIsInt
   res= cast(round(filtered),imType);
else
   res=cast(filtered,imType);
end





