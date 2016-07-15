function srInMeshOut = filterSrBact(srInMesh,varargin)

parInfo = srInMesh{1}.localizations.parInfo;
photCol = findSRField(parInfo,'semantic','emission strength');
bgCol   = findSRField(parInfo,'semantic','fit residue chi square value');
chi2Col = findSRField(parInfo,'semantic','local background');
zCol = findSRField(parInfo,'semantic','position in sample space in z dimension');
frameCol =  findSRField(parInfo,'semantic','frame number');

filterPhot = false;
filterBg = false;
filterChi2 = false;
filterZ = false;
filterFrame = false;
plotOn = false;

ii =1;
while ii <= numel(varargin)
   if strcmp(varargin{ii},'FilterPhot')
      filterPhot = true;
      photLim= varargin{ii+1};
      ii = ii +2;
   elseif strcmp(varargin{ii},'FilterBg')
      filterBg = true;
      bgLim= varargin{ii+1};
      ii = ii +2;
   elseif strcmp(varargin{ii},'FilterChi2')
      filterChi2 = true;
      chi2Lim= varargin{ii+1};
      ii = ii +2;
   elseif strcmp(varargin{ii},'FilterZ')
      filterZ = true;
      zLim= varargin{ii+1};
      ii = ii +2;
   elseif strcmp(varargin{ii},'PlotOn')
      plotOn = true;
      nBin = varargin{ii+1};
      ii = ii + 2;
   elseif strcmp(varargin{ii},'FilterFrame')
     filterFrame = true;
     frameLim = varargin{ii+1};
      ii = ii + 2;
   else 
      ii = ii+1;
   end
end

if plotOn;
   nCell = numel(srInMesh);
   phot = [];
   bg = [];
   z=[];
   chi2 = [];
   for ii = 1:nCell
      phot =  [phot; srInMesh{ii}.localizations.data(:,photCol)];
      bg =    [bg;srInMesh{ii}.localizations.data(:,bgCol)];
      chi2 =  [chi2;srInMesh{ii}.localizations.data(:,chi2Col)];
      z =  [z;srInMesh{ii}.localizations.data(:,zCol)];
   end

   subplot(4,1,1);
   hist(phot,nBin);
   xlabel('Photon count');
   subplot(4,1,2);
   hist(bg,nBin);
   xlabel('Local background');
   subplot(4,1,3);
   hist(chi2,nBin);
   xlabel('Chi2');
   subplot(4,1,4);
   hist(z,nBin);
   xlabel('Z');
end

if filterPhot | filterBg | filterChi2 | filterZ | filterFrame
   nCell = numel(srInMesh);
   for ii = 1:nCell
      phot =  [srInMesh{ii}.localizations.data(:,photCol)];
      bg =    [srInMesh{ii}.localizations.data(:,bgCol)];
      chi2 =  [srInMesh{ii}.localizations.data(:,chi2Col)];
      z =  [srInMesh{ii}.localizations.data(:,zCol)];
      frame =  [srInMesh{ii}.localizations.data(:,frameCol)];
      cellOk = ones(size(phot));
      if filterPhot 
         cellOk = cellOk & (phot> photLim(1) & phot< photLim(2));
      end 
      if filterBg
         cellOk = cellOk & (bg > bgLim(1) & bg < bgLim(2));
      end 
      if filterChi2
         cellOk = cellOk & (chi2> chi2Lim(1) & chi2< chi2Lim(2));
      end 
      if filterZ
         cellOk = cellOk & (z> zLim(1) & z< zLim(2));
      end 
      if filterFrame
         cellOk = cellOk & (frame> frameLim(1) & frame< frameLim(2));
      end 

      srInMesh{ii}.l(~cellOk) = [];
      srInMesh{ii}.d(~cellOk) = [];
      srInMesh{ii}.position(~cellOk) = [];
      srInMesh{ii}.localizations.data(~cellOk(:),:) = [];
   end
end

srInMeshOut = srInMesh;
