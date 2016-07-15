function res = valleySegmentation2(im,threshminlevel,opennum,erodenum,edgeSigmaV,valleythresh1,valleythresh2,invertimage,thresFactorM)
  %microbetrackers valley segmentation routine, for 
  %combination with marker based segmentation
  %p.threshminlevel

  if invertimage
     im = max(im(:))-im;
  end

  imge = img2imge(im,erodenum);
  imge16 = img2imge16(im,erodenum);

  thres = graythreshreg(imge,threshminlevel);
  edgemode = 'valley';
  edgeSigmaL=[];
  res= getRegions(edgemode,imge,thres,thresFactorM,opennum,imge16,edgeSigmaL,edgeSigmaV,valleythresh1,valleythresh2);
  %res = bwlabel(res>0,4);
end
%-------------------------------------------------------------------
%-------------------------------------------------------------------
function res = graythreshreg(img,varargin)
  % threshold calculated in a regionSelectionRect region
  if ~isempty(varargin), flevel=varargin{1}; else flevel=0; end
  %global regionSelectionRect
  regionSelectionRect=[];
  sz = size(img);
  if isempty(regionSelectionRect)
      res = graythresh2(img(ceil(sz(1)*0.05):floor(sz(1)*0.95),ceil(sz(2)*0.05):floor(sz(2)*0.95),1));
  else
      res = graythresh2(imcrop(img,regionSelectionRect));
  end
  function b=graythresh2(a)
      if flevel>0
          c = reshape(a,1,[]);
          c = sort(c);
          level = c(ceil(min(flevel,1)*length(c)));
          b = graythresh(c(c>=level));
      else
          b = graythresh(a);
      end
  end
end

%-------------------------------------------------------------------
%-------------------------------------------------------------------
function res=getRegions(edgemode,im,thres,tfactor,opennum,varargin)
% This function determines and labels regions on an image (im) using edge
% detection and thresholding (tractor - factor compared to automatic
% threshols, larger values - smaller regions)
% 
% mode 'none' - no edge/valley detection, thresholding only (params: im & tfactor)
% mode 'log' - thresholding + LoG edge detection (params: edgeSigmaL)
% mode 'valley' - thresholding + valley detection (params: edgeSigmaV,valleythresh1,valleythresh2)
% mode 'logvalley' - thresholding + LoG + valley (params: all)

  if ~isempty(varargin)
      im16 = varargin{1};
      edgeSigmaL = varargin{2};
  end
  if length(varargin)>=3
      edgeSigmaV = varargin{3};
      valleythresh1 = varargin{4};
      valleythresh2 = varargin{5};
  end

  thres2 = max(0,min(1,thres*tfactor));
  if thres*tfactor>=1, gdisp('Warning: threshold exceeded 1, no thresholding performed. thresFactorM may be too high.'); end
  mask = im2bw(im,thres2);

  edg = logvalley(im16,edgemode,edgeSigmaL,edgeSigmaV,valleythresh1,valleythresh2);
  edg = bwmorph(edg,'clean');
  edg = bwmorph(edg,'bridge');
  imgProc = (1-edg).*mask;

  seo = strel('disk',max(floor(opennum),0));
  imgProc = imopen(imgProc,seo);
  imgProc(:,[1 end])=0;
  imgProc([1 end],:)=0;
  res = bwlabel(imgProc,4);
end
%-------------------------------------------------------------------
%-------------------------------------------------------------------
function g = logvalley(img,mode,sigmaL,sigmaV,thresh1,thresh2)
  % a combination of the 'valley detection' and 'LoG' algorithms
  % 
  % mode = 'none', 'log', 'valley', 'logvalley' or 0,1,2,3
  % thresh2 is a 'hard' threshold, thresh1 - a 'soft' one, meaning that the
  % pixel is only detected as a part of a valley if it's adjasent to a
  % pixel above the hard threshold.

  if isnumeric(mode)
      if mode==0, mode='none';
      elseif mode==1, mode='log';
      elseif mode==2, mode='valley';
      elseif mode==3, mode='logvalley';
      end
  end
  img = 1-im2single(img);
  [m,n] = size(img);

  if strcmp(mode,'log') || strcmp(mode,'logvalley')
      % LoG edge  detection
      if sigmaL<0.01, sigmaL=0.01; end
      fsize = ceil(sigmaL*3)*2 + 1;  % choose an odd fsize > 6*sigmaL;
      op = fspecial('log',fsize,sigmaL); 
      op = op - sum(op(:))/numel(op); % make the op to sum to zero
      b = imfilter(img,op,'replicate');
      e = b>0;
      se = strel('arbitrary',ones(3));
      e = imdilate(e&~bwmorph(e,'endpoints'),se)&e;
      if isempty(thresh1), thresh1 = .75*mean2(abs(b)); end
  end
  if strcmp(mode,'valley') || strcmp(mode,'logvalley') || strcmp(mode,'flogvalley')
      % Valley detection
      se = strel('arbitrary',ones(3));
      rr = 2:m-1; 
      cc = 2:n-1;
      if isempty(thresh1), thresh1 = 0; end
      if length(thresh1)>1, thresh1 = thresh1(rr,cc); end
      if isempty(thresh2), thresh2 = 0; end
      thresh2 = max(thresh2,0);
      thresh1 = min(max(thresh1,0),thresh2);
      if sigmaV>0
          fsize = ceil(sigmaV*3)*2 + 1;
          op = fspecial('gaussian',fsize,sigmaV); 
          op = op/sum(sum(op));
          img = imfilter(img,op,'replicate');
      end
      f = zeros(m,n);
      f(rr,cc) = max(max(min(img(rr,cc-1)-img(rr,cc),img(rr,cc+1)-img(rr,cc)), ...
                         min(img(rr-1,cc)-img(rr,cc),img(rr+1,cc)-img(rr,cc))), ...
                     max(min(img(rr-1,cc-1)-img(rr,cc),img(rr+1,cc+1)-img(rr,cc)), ...
                         min(img(rr+1,cc-1)-img(rr,cc),img(rr-1,cc+1)-img(rr,cc))));
      fmean = mean(f(f>quantile(f(f>0),0.99)));
      f1 = f>fmean*thresh1;
      f2 = f>fmean*thresh2;
      for i=1:4, f2 = f1 & imdilate(f2,se); end
  end
  if strcmp(mode,'log') || strcmp(mode,'flog')
      g = e;
  elseif strcmp(mode,'valley')
      g = f2;
  elseif strcmp(mode,'logvalley') || strcmp(mode,'flogvalley')
      g = e | f2;
  else
      g = false(m,n);
  end
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

