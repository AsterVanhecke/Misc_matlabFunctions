function [rCM,rMax, intRadial] = getRadialProfile(im,xCM,yCM,maxRad,varargin)

plotOn = false;
ii = 1;
while ii <=numel(varargin)
   if strcmp(varargin{ii},'PlotOn')
      plotOn = true;
      ii = ii+1;
   else
      ii = ii + 1;
   end
end

intRadial = zeros(maxRad,1);

if ~isnan(xCM) && ~isnan(yCM)
   [sizey,sizex] = size(im);

   distToEdgeX = min(xCM,sizex-xCM);
   distToEdgeY = min(yCM,sizey-yCM);
   rMaxObs = floor(min(distToEdgeX,distToEdgeY));

   [X,Y] = meshgrid(1:sizex,1:sizey);

   X = X - xCM;
   Y = Y - yCM;

   [THETA,RHO] = cart2pol(X,Y);


   for ii = 1:rMaxObs
      intRadial(ii) =sum(sum(im(RHO >= ii-.5 & RHO <= ii+.5 )));
   end

   r = [1:maxRad]';
   rCM = sum(r.*intRadial)/ sum(intRadial);
   [iMax, idx ]= max(intRadial);
   rMax = r(idx);

   if plotOn
      r = 1:maxRad;
      plot(r,intRadial);
      xlabel('Radius, pix');
      ylabel('Radial intensity');
   end
else
   rCM = NaN;
   rMax = NaN;
end
