function m = mTrack()
  % Ripped out of the microbe Tracker GUI
  %
  %Available functions:
  %loadimages(channel,folder) -- loads a series of TIFF images from <folder> to <channel>, which is 1 for phase, 3 for signal 1, 4 for signal 2
  %loadimagestack(channel,file) -- loads a stack of images from <file> to <channel> (same channel numbering)
  %loadparam(file) -- load a set of saved parameters from <file>
  %param = loadmesh(filename) -- loads mesh from <filename> (<filename> should be a full path or a name in the current folder)
  %parseparams(param) -- initializes the parameter set with the string <param> from loadmesh command
  %alignfrm -- aligns the frames loaded at the moment
  %savealign(filename) -- saves the alignment data to file <filename>
  %loadalign(filename) -- loads the alignment data from file <filename>
  %subtractbgr(channels,range,invert) -- subtract background from <channels> = <ch1> or <channels> = [<ch1> <ch2>], 3 for signal 1, 4 for signal 2
  %process(range,mode,lst,addsig,addas,savefile,fsave,saveselect,region,shiftfluo) -- main processing function,
  %  range - range of frames (an array of two numbers or a single number) on which the processing will be made or an empty array for all frames, e.g. [1 10] for frames from 1 to 10, [] for all frames.
  %  mode - 1-Timelapse, 2-1st independent frame, 3-Independent frames, 4-Reuse meshes.
  %  lst - list of cells to be processed, defined on the frame preceeding the first frame in the range (if detection is performed) or on any frame (in Reuse meshes mode). Use the empty array [] to process all cells.
  %  addsig - an array of four numbers 0 or 1, indicating whether to add signal to each channel (as above: 1-phase, 3-signal 1, 4-signal 2).
  %  addas - a cell array of values indicating the names under which the signals will be saved. Default {}, which is equivalent to {0,-,1,2}. Each numeric value X is converted to signalX, each string value remains unchanged. I.e. {0,0,'s1',200} will save the profile of the phase signal to field "signal0", signal 1 to field "s1", signal 2 to field "signal200".
  %  savefile - filename to save the analysis while processing (such as 'tempresults.mat' in the main window). Provide an enpty string '' if no saving is necessary.
  %  fsave - "frequency" of saveing, i.e. the number of frames per one saving.
  %  saveselect - 1 to save only selected cells and 0 to save all cells.
  %  region - region to process (an array of four numbers, left-top- width-height), [] for the whole image.
  %savemesh(filename,list,mode,range) -- save the mesh, <list> - list of selected cells, <mode> = true - save selected, false - save all, <range> - frame range
  %assignSegment(x,y,mesh,box,coordOffset) - Assign segment number in which spot is located inside cell
  %[XIn] = dblTformInv(tformDbl,XBase)
  %dblTform_InToBase = getDblTform(In,Base) 

  % initalise some global variables as per microbeTracker gui
  initMTrack;

  %add the functions to the data structure
  m.loadimages	    =  @loadimages;			
  m.loadimagestack  =  @loadimagestack;			
  m.loadparam	    =  @loadparam;			
  m.loadmesh	    =  @loadmesh;			
  m.parseparams	    =  @parseparams;			
  m.alignfrm	    =  @alignfrm;			
  m.savealign	    =  @savealign;			
  m.loadalign	    =  @loadalign;			
  m.subtractbgr	    =  @subtractbgr;			
  m.process	    =  @process;			
  m.savemesh	    =  @savemesh;			
  m.assignSegment   =  @assignSegment;
  m.dblTformInv     = @dblTformInv;
  m.getDblTform     = @getDblTform;
  m.returnparam = @returnparam;
end
			
function initModel			
  % used in "process" function to define a set of array (globals
  % variables) according to the chosen algorithm and some other
  % parameters. These arrays are used later in the corresponding
  % PDM-based model (algoritms 2 and 3)
  global coefPCA coefS mCell N weights dMax p
  if p.algorithm==2
    load(p.trainingFile,'coefPCA','latPCA','mCell')
    N = size(coefPCA,1)/2;
    coef4 = [mCell(:,1);mCell(:,2)];
    coef4max = sqrt(sum(coef4.^2));
    coef4 = coef4/coef4max;
    coefPCA = [[ones(N,1) zeros(N,1);zeros(N,1) ones(N,1)]/sqrt(N) coef4 coefPCA(:,1:p.Nkeep)];
    weights = 1*[1;1;1;sqrt(latPCA(1:p.Nkeep))/sqrt(latPCA(1))];
    dMax = (1/1)*sqrt(max((mCell(:,1)-mean(mCell(:,1))).^2+(mCell(:,2)-mean(mCell(:,2))).^2));
  elseif p.algorithm==3
    % scaleFactor
    load(p.trainingFile,'coefPCA','latPCA','mCell')
    if isfield(p,'scaleFactor'), mCell=mCell.*p.scaleFactor;
    else p.scaleFactor = 1;
    end
    N = size(coefPCA,1)/2;
    coef4 = sin([mCell(:,1);0*mCell(:,2)]/max(mCell(:,1))*pi/2);
    coef4max = sqrt(sum(coef4.^2));
    coef4 = coef4/coef4max;
    coefS = [ones(N,1) zeros(N,1);zeros(N,1) ones(N,1)]/sqrt(N);
    coefPCA = [coef4 coefPCA(:,1:p.Nkeep)];
    weights = 1*[1;1;1;0.98;0.98;0.98*sqrt(latPCA(1:p.Nkeep))/sqrt(latPCA(1))/2]; %08/18/08
    dMax = (1/1)*sqrt(max((mCell(:,1)-mean(mCell(:,1))).^2+(mCell(:,2)-mean(mCell(:,2))).^2));
  end
end

%% "process" function
function process(range,mode,lst,addsig,addas,savefile,fsave,saveselect,processregion,shiftfluo)
  % Main processing function, decides on which frames and what type of
  % cell detection or cell analysis to do, to be called by "Detection &
  % analysis" buttons or from eacher variant of the batch mode
  % 
  % "range" - range of frames to run, can be []-all frames
  % "mode" - 1-tlapse, 2-1st ind, 3-all ind, 4-reuse
  % "lst" - list of cells on the frame previous to range(1)
  % "addsig" - [0 0 0 0]=[0 0 0]=[0 0]=0=[]-no signal, 1st one phase, etc. 
  % "addas" - default {}={0,-,1,2}, if numeric X, creates signalX, else - X
  % "savefile" - filename to save to
  % "fsave" - frequance of saving, n-once per n frames, 0-never, []-end
  
  global p cellList imsizes
  
  if isempty(p) || ~isfield(p,'algorithm'), gdisp('Error: parameters not initialized.'); return; end
  
  initModel
  
  if mode==4 && ~isempty(range) && (length(cellList)<range(1)), gdisp('Processing error: selected frame out of range'); return; end
  if mode==4 && isempty(range) && isempty(cellList{1}), gdisp('Processing error: for all frames, the 1st frame must not be empty'); return; end
  if mode~=4 && ~isempty(lst) && (isempty(range) || range(1)==1), gdisp('Processing error: selected cells regime does not work on the 1st frame'); return; end
  
  if ~isempty(fsave) && (fsave==0 || isempty(savefile)), savemode = 1; % never save
  elseif isempty(fsave), savemode = 2; % save at the end
  else savemode = 3; % save on some steps
  end
  if length(range)==1, range = [range range]; end
  if isempty(range), range = [1 imsizes(end,3)]; end
  
  time1 = clock;
  listtorun = lst;
  for frame=range(1):range(2)
    if isempty(lst) && mode~=4, cellList{frame}={}; end
    if mode==3 || (mode~=4 && frame==1) || (mode==2 && frame==range(1)) || isempty(cellList) ||...
        (frame>=2 &&(length(cellList)<frame-1||(isempty(cellList{frame-1})))&&mode~=4)
      if p.algorithm==1
        processFrameOld(frame,[],false,processregion);
      elseif ismember(p.algorithm,[2 3 4])
        processFrameI(frame,ismember(mode,[2 3]),processregion);
        joincells(frame,[]); 
      end
    elseif mode==1 || mode==2
      shiftmeshes(frame,1)
      if isempty(lst), listtorun=[]; end
      if p.algorithm==1, listtorun = processFrameOld(frame,listtorun,true,[]); end
      if ismember(p.algorithm,[2 3 4]), listtorun = processFrame(frame,listtorun); end
      shiftmeshes(frame,-1)
    elseif mode==4
      if isfield(p,'joinWhenReuse') && p.joinWhenReuse, 
        joincells(frame,[]); 
      end
      if frame~=range(1)
        listtorun = selNewFrame(listtorun,frame-1,frame);
      end
    end
    if length(cellList)>=frame
      for i=1:length(cellList{frame})
        if ~isempty(cellList{frame}{i})
          cellList{frame}{i}=getextradata(cellList{frame}{i});
        end
      end
    end
    if isempty(lst), listtorun=[]; end
    for i=1:length(addsig)
      if addsig(i)
        if length(addas)<i, addas0 = max(i-2,0); else addas0 = addas{i}; end
        addSignalToFrame(frame,i,addas0,listtorun,p.sgnResize,p.approxSignal,shiftfluo);
      end
    end
    if savemode==3 && mod(frame,fsave)==0, savemesh(savefile,selNewFrame(lst,frame,1),saveselect,[]); end
    time2 = clock;
    gdisp(['frame ' num2str(frame) ' finished, elapsed time ' num2str(etime(time2,time1)) ' s']);
    gdisp(' ');
  end
  if savemode==2 || (savemode==3 && mod(frame,fsave)~=0), savemesh(savefile,lst,saveselect,[]); end
end

% This function has been taken out from the "process" function to save space
function shiftmeshes(frame,dir)
  % shifts the position of every cell on the indicated frame in "cellList"
  % structure according to the amount indicated in the "shiftframes" array
  % (obtained from image alignment) in the indicated direction "dir" (1 -
  % direct, -1 - backward). This function is used to prepare cell
  % positions from the previous frame to be used as the initial guess in
  % an active coutour routine.
  global p cellList shiftframes imsizes
  if isempty(who('shiftframes')) || isempty(shiftframes) || length(shiftframes.x)~=imsizes(1,3) || ~ismember(dir,[-1 1]) , return; end
  for c=1:length(cellList{frame-1})
    if ~isempty(cellList{frame-1}{c}) && isfield(cellList{frame-1}{c},'model') && ismember(p.algorithm,2)
      cellList{frame-1}{c}.model(2) = cellList{frame-1}{c}.model(2)+dir*shiftframes.y(frame-1)-dir*shiftframes.y(frame);
      cellList{frame-1}{c}.model(3) = cellList{frame-1}{c}.model(3)+dir*shiftframes.x(frame-1)-dir*shiftframes.x(frame);
    elseif ~isempty(cellList{frame-1}{c}) && isfield(cellList{frame-1}{c},'model') && ismember(p.algorithm,3)
      cellList{frame-1}{c}.model(1) = cellList{frame-1}{c}.model(1)+dir*shiftframes.y(frame-1)-dir*shiftframes.y(frame);
      cellList{frame-1}{c}.model(2) = cellList{frame-1}{c}.model(2)+dir*shiftframes.x(frame-1)-dir*shiftframes.x(frame);
    elseif ~isempty(cellList{frame-1}{c}) && isfield(cellList{frame-1}{c},'model') && ~isempty(cellList{frame-1}{c}.model) && ismember(p.algorithm,4)
      cellList{frame-1}{c}.model(:,1) = cellList{frame-1}{c}.model(:,1)+dir*shiftframes.y(frame-1)-dir*shiftframes.y(frame);
      cellList{frame-1}{c}.model(:,2) = cellList{frame-1}{c}.model(:,2)+dir*shiftframes.x(frame-1)-dir*shiftframes.x(frame);
    elseif ~isempty(cellList{frame-1}{c}) && p.algorithm==1 && isfield(cellList{frame-1}{c},'mesh') && length(cellList{frame-1}{c}.mesh)>=4
      cellList{frame-1}{c}.mesh(:,[1 3]) = cellList{frame-1}{c}.mesh(:,[1 3])+dir*shiftframes.y(frame-1)-dir*shiftframes.y(frame);
      cellList{frame-1}{c}.mesh(:,[2 4]) = cellList{frame-1}{c}.mesh(:,[2 4])+dir*shiftframes.x(frame-1)-dir*shiftframes.x(frame);
    elseif ~isempty(cellList{frame-1}{c}) && p.algorithm==1 && isfield(cellList{frame-1}{c},'contour')
      cellList{frame-1}{c}.contour(:,1) = cellList{frame-1}{c}.contour(:,1)+dir*shiftframes.y(frame-1)-dir*shiftframes.y(frame);
      cellList{frame-1}{c}.contour(:,2) = cellList{frame-1}{c}.contour(:,2)+dir*shiftframes.x(frame-1)-dir*shiftframes.x(frame);
    end
  end
end

%% Cell statistics and additional data

function str=getextradata(str)
  % calculates geometrical properties of detected meshes in structure
  % "str", corresponding to a single cell in "cellList", recording to the
  % same structure. The properties are: "steplength", "steparea",
  % "stepvolume" (length, area, and volume of each segment), "length",
  % "area", "volume" (foe the whole cell), "lengthvector" - coordinates
  % along cell centerline of the centers of each segment.
  if isempty(str), return; end
  if isfield(str,'mesh') && length(str.mesh)>1
    mesh = str.mesh;
    lng = size(mesh,1)-1;
    if ~isfield(str,'polarity'), str.polarity=0; end

    %length
    str.steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
      mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
    str.length = sum(str.steplength);

    str.lengthvector = cumsum(str.steplength)-str.steplength/2;

    % area
    str.steparea = [];
    for i=1:lng, str.steparea=[str.steparea;polyarea([mesh(i:i+1,1);mesh(i+1:-1:i,3)],[mesh(i:i+1,2);mesh(i+1:-1:i,4)])]; end
    str.area = sum(str.steparea);

    % volume
    d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
    str.stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*str.steplength*pi/4;
    str.volume = sum(str.stepvolume);
  elseif isfield(str,'contour') && length(str.contour)>1
    contour = str.contour;
    lng = size(contour,1);
    str.length = sqrt(max(max( (repmat(contour(:,1),1,lng)-repmat(contour(:,1)',lng,1)).^2 + ...
                   (repmat(contour(:,2),1,lng)-repmat(contour(:,2)',lng,1)).^2)));
    str.area = polyarea(contour(:,1),contour(:,2));
  end
end

function d=edist(x1,y1,x2,y2)
  % complementary for "getextradata", computes the length between 2 points
  d=sqrt((x2-x1).^2+(y2-y1).^2);
end

%% Cell orientations

function lst = removeorientationall(lst,cell)
  % removes orientation for the given cell on each frame, updating the
  % list "lst" (clone of "cellList")  
  for frame=1:length(lst)
    if cell<=length(lst{frame}) && ~isempty(lst{frame}{cell})
      lst{frame}{cell}.polarity = 0;
    end
  end
end

function lst = reorientall(lst,cell,flag)
  % reorients the given cell ("cell") on each frame in structure lst if
  % flag==true, but setting the polarity porameter to true in either case
  % 
  tflag = false;
  for frame=1:length(lst)
    if cell<=length(lst{frame}) && ~isempty(lst{frame}{cell})
      % terminate if the cell just divided (except on the 1st frame)
      if tflag && lst{frame}{cell}.birthframe==frame
        return
      end
      tflag = true;
      % reorient
      if flag
        lst{frame}{cell} = reorient(lst{frame}{cell});
      end
      % set the polarity variable
      lst{frame}{cell}.polarity = 1;
    end
  end
end

function str2=reorient(str)
  % this function reorients the data corresponding to one cell on a 
  % single frame (without changing the polarity variable)
  % 
  names = fieldnames(str);
  % flip the orientation of the data arrays (of every vertical array)
  for i=1:length(names)
    t=['str2.' names{i} '=flipud(str.' names{i} ');'];
    eval(t);
  end
  % if mesh is flipped, preserve the clockwise orientation of the outline
  if ismember('mesh',names) && numel(str2.mesh)>1
    str2.mesh = str2.mesh(:,[3 4 1 2]);
  end
  % change the models for particular algorithms
  if ismember('model',names) && ismember('algorithm',names)
    if ismember(str.algorithm,[1 2])
      str2.model(1) = mod(str.model(1) + pi,2*pi);
    elseif ismember(str.algorithm,3)
      str2.model(3) = mod(str.model(3) + pi,2*pi);
    elseif ismember(str.algorithm,4)
      s = size(str.model,1);
      s2 = floor(s/2);
      if s>1
        str2.model = [str.model(s2+1:s,:);str.model(1:s2,:)];
      end
    end
  end
  % reorient spots data (assuming the spots structure name starts with
  % 'spots' and all spots sub-arrays are horizontal)
  for i=1:length(names)
    if length(names{i})>=5 && strcmp('spots',names{i}(1:5))
      names2 = {};
      eval(['names2 = fieldnames(str.' names{i} ');'])
      if ismember('positions',names2)
        eval(['str2.' names{i} '.l = str.length-str.' names{i} '.l;'])
        eval(['str2.' names{i} '.positions = size(str.mesh,1)+1-str.' names{i} '.positions;'])
        eval(['names2 = fieldnames(str.' names{i} ');'])
        for j=1:length(names2)
          t=['str2.' names{i} '.' names2{j} '=fliplr(str2.' names{i} '.' names2{j} ');'];
          eval(t);
        end
      end
    end
  end
  % reorient lengthvector (the only aray for which flip is insufficient)
  if ismember('lengthvector',names) && ismember('length',names)
    str2.lengthvector = str2.length-str2.lengthvector;
  end
end

%% Aligning global functions

function loadalign(filename)
  global shiftframes
  % This function loads aligning data from file <filename>
  % The file must exist and be a .mat file
  % Intended use: with "alignphaseframes" callback & batch files
  loaded = load(filename,'shiftframes');
  if ~isfield(loaded,'shiftframes'), gdisp('This file does not contain alignment data'); return; end
  shiftframestmp = loaded.shiftframes;
  if ~isstruct(shiftframestmp), gdisp('This file does not contain alignment data'); return; end
  if ~isfield(shiftframestmp,'x') || ~isfield(shiftframestmp,'y'), gdisp('This file does not contain alignment data'); return; end
  shiftframes.x = shiftframestmp.x;
  shiftframes.y = shiftframestmp.y;
  gdisp('Alignment data loaded')
end

function savealign(filename)
  global shiftframes %#ok<NUSED>
  % This function saves aligning data to file <filename>
  % The file must be a .mat file
  % Intended use: with "alignphaseframes" callback & batch files
  save(filename,'shiftframes');
  gdisp('Alignment data saved')
end

function alignfrm
  % This function alignes phase images
  % The result is stored in "shiftframes" structure with x & y fields of the main function 
  % Intended use: with "alignphaseframes" callback & batch files
  global rawPhaseData p shiftframes
  
  if checkparam(p,'aligndepth'), gdisp('Images not aligned: parameter "aligndepth" not provided.'); return; end
  if isempty(rawPhaseData), gdisp('Images not aligned: no phase images loaded.'); return; end
  [shiftframes.x,shiftframes.y]=alignframes(rawPhaseData,p.aligndepth);
  gdisp('Images aligned')
end

%% Initialization global functions

function cleardata
  global rawPhaseData rawS1Data rawS2Data cellList cellListN selectedList shiftframes imsizes handles %#ok<NUSED>
  % close open windows
  if exist('handles','var') && isstruct(handles)
    fields = fieldnames(handles);
    for i=1:length(fields)
      eval(['cfield = handles.' fields{i} ';'])
      if ishandle(cfield)
        delete(cfield);
      elseif isstruct(cfield)
        fields2 = fieldnames(cfield);
        for k=1:length(cfield)
          for j=1:length(fields2)
            eval(['if ishandle(cfield(' num2str(k) ').' fields2{j} '), delete(cfield(' num2str(k) ').' fields2{j} '); end'])
          end
        end
      end
    end
  end
  % cleas variables
  clear global rawPhaseData rawS1Data rawS2Data cellList cellListN selectedList shiftframes imsizes handles
end

%% Batch global functions

function batchtextrun_glb(s)
  global rawPhaseData rawS1Data rawS2Data cellList cellListN p
  eval(s);
end

%% Images global function

function subtractbgr(channels,range,varargin)
  % background subrtaction routine:
  % "channels" - list of channels (3-signal 1, 4-signal 2)
  % "range" - [first_frame last_frame] (empty = all frames)
  % "invert" (optional) - invert the "phase" image
  
  global rawPhaseData rawS1Data rawS2Data imsizes se p
  
  if checkparam(p,'bgrErodeNum') || ~strcmp(class(p.bgrErodeNum),'double')
    gdisp('Background subtraction failed: parameter "bgrErodeNum" not provided.');
    return
  end
  if length(varargin)>=1
    invert = varargin{1};
  else
    invert = false;
  end
  if isempty(channels), channels = [3 4]; end
  if min(imsizes(channels,1))<2, return; end
  if isempty(rawPhaseData), gdisp('Background subtraction failed: no phase images loaded'); return; end
  if isempty(range), range = [1 10000]; end
  if length(range)==1, range = [range range]; end
  f = channels*0;
  for g=1:length(channels)
    if channels(g)==3
      if imsizes(3,1)<2, continue;
      else crange=[max(1,range(1)) min(imsizes(3,3),range(2))];
      end
    elseif channels(g)==4
      if imsizes(4,1)<2, continue; % isempty(who('rawS2Data')) || isempty(rawS2Data), 
      else crange=[max(1,range(1)) min(imsizes(4,3),range(2))];
      end
    else
      continue
    end
    if imsizes(1,3)==1 && crange(2)>crange(1)
      gdisp('Background subtraction warning: using single phase contrast image for multiple signal images');
    elseif imsizes(1,3)<crange(2)
      gdisp('Background subtraction error: the number of phase contrast and signal images not matching');
      continue
    end
    for i=crange(1):crange(2)
      f(g) = 1;
      if imsizes(1,3)==1
        imgP = rawPhaseData(:,:,1);
      else
        imgP = rawPhaseData(:,:,i);
      end
      if isempty(imgP), continue; end
      if channels(g)==3, img = rawS1Data(:,:,i); end
      if channels(g)==4, img = rawS2Data(:,:,i); end
      if size(img,1)~=size(imgP,1) || size(img,2)~=size(imgP,2)
        gdisp('Background subtraction error: the images are of different size.')
        break
      end
      if invert, imgP = max(max(imgP))-imgP; end
      thres = graythreshreg(imgP,p.threshminlevel);
      mask = im2bw(imgP,thres);
      for k=1:p.bgrErodeNum, mask = imerode(mask,se); end
      bgr = mean(img(mask));
      img0 = int32(img);
      img1 = img0-bgr;
      img2 = max(0,img1);
      img = uint16(img2);
      if channels(g)==3, rawS1Data(:,:,i) = img; end
      if channels(g)==4, rawS2Data(:,:,i) = img; end
      if mod(i,5)==0, gdisp(['Subtracting backgroung from signal ' num2str(channels(g)-2) ', frame ' num2str(i)]); end
    end
  end
  if sum(f)>1
    gdisp(['Subtracting backgroung completed from ' num2str(sum(f)) ' channels']);
  elseif sum(f)==1
    gdisp(['Subtracting backgroung completed from signal ' num2str(channels(f)-2)]);
  end
end

function bgr = phasebgr(img,thres)
  global se p
  mask = ~im2bw(img,thres);
  for k=1:p.bgrErodeNum, mask = imerode(mask,se); end
  bgr = mean(img(mask));
end

function res = loadimagestack(n,filename)
  % loads image stacks using Bioformats
  % "n" - channel (1-phase, 2-extra, 3-signal1, 4-signal2)
  % "filename" - stack file name
  % 
  % This function differs from the same named function in the suite
  % folder, which is used by SpotFinder tools. The differences are: (1)
  % loading files directly into the global structures to decrease memory
  % requirements; (2) use of a waitbar. The same function as here is also
  % used by aligntool.
  % 
  global rawPhaseData rawS1Data rawS2Data imsizes imageLimits
  bformats = checkbformats(0);
  res = false;
  if n==1
    str='rawPhaseData';
  elseif n==3
    str='rawS1Data';
  elseif n==4
    str='rawS2Data';
  else
    gdisp('Error loading images: channel not supported');
  end
  if (length(filename)>4 && strcmpi(filename(end-3:end),'.tif')) || (length(filename)>5 && strcmpi(filename(end-4:end),'.tiff'))
    % loading TIFF files
    try
      info = imfinfo(filename);
      numImages = numel(info);
      w = waitbar(0, 'Loading images, please wait...');
      lng = info(1).BitDepth;
      if lng==8
        cls='uint8';
      elseif lng==16
        cls='uint16';
      elseif lng==32
        cls='uint32';
      else
        gdisp('Error in image bitdepth loading multipage TIFF images: no images loaded');return;
      end
      eval([str '=zeros(' num2str(info(1).Height) ',' num2str(info(1).Width) ',' num2str(numImages) ',''' cls ''');'])
      for i = 1:numImages
        img = imread(filename,i,'Info',info);
        eval([str '(:,:,' num2str(i) ')=img;'])
        waitbar(i/numImages, w);
      end
      close(w)
      eval(['imageLimits{n} = 2^' num2str(lng) '*mean(stretchlim(' str ',[0.0001 0.9999]),2);']);
      eval(['imsizes(n,:) = [size(' str ',1) size(' str ',2) size(' str ',3)];']);
      gdisp(['Loaded ' num2str(imsizes(n,3)) ' images from a multipage TIFF'])
      updateimsizes
      res = true;
    catch
      gdisp('Error loading multipage TIFF images: no images loaded');
    end
  elseif bformats
    % loading all formats other than TIFF
    try
      breader = loci.formats.ChannelFiller();
      breader = loci.formats.ChannelSeparator(breader);
      breader = loci.formats.gui.BufferedImageReader(breader);
      breader.setId(filename);
      numSeries = breader.getSeriesCount();
      if numSeries~=1, gdisp('Incorrect image stack format: no images loaded'); return; end; 
      breader.setSeries(0);
      wd = breader.getSizeX();
      hi = breader.getSizeY();
      shape = [wd hi];
      numImages = breader.getImageCount();
      if numImages<1, gdisp('Incorrect image stack format: no images loaded'); return; end;
      nBytes = loci.formats.FormatTools.getBytesPerPixel(breader.getPixelType());
      if nBytes==1
        cls = 'uint8';
      else
        cls = 'uint16';
      end
      eval([str '=zeros(' num2str(hi) ',' num2str(wd) ',' num2str(numImages) ',''' cls ''');'])
      w = waitbar(0, 'Loading images, please wait...');
      for i = 1:numImages
        img = breader.openImage(i-1);
        pix = img.getData.getPixels(0, 0, wd, hi, []);
        arr = reshape(pix, shape)';
        if nBytes==1
          arr2 = uint8(arr/256);
        else
          arr2 = uint16(arr);
        end
        eval([str '(:,:,' num2str(i) ')=arr2;'])
        waitbar(i/numImages, w);
      end
      close(w)
      % eval(['cls = class(' str ');']);
      if strcmp(cls,'uint8'), lng=8; elseif strcmp(cls,'uint16'), lng=16; elseif strcmp(cls,'uint32'), lng=32;
        else gdisp('Error in image bitdepth loading images using BioFormats: no images loaded');return; end
      eval(['imageLimits{n} = 2^' num2str(lng) '*mean(stretchlim(' str ',[0.0001 0.9999]),2);']);
      eval(['imsizes(n,:) = [size(' str ',1) size(' str ',2) size(' str ',3)];']);
      gdisp(['Loaded ' num2str(imsizes(n,3)) ' images using BioFormats'])
      updateimsizes
      res = true;
    catch
      gdisp('Error loading images using BioFormats: no images loaded');
    end
  else % unsupported format of images
    gdisp('Error loading images: the stack must be in TIFF format for BioFormats must be loaded');
  end
end

function loadimages(n,folder)
  % loads TIFF (image) files into the specified channel:
  % "n" - channel (1-phase, 2-extra, 3-signal1, 4-signal2)
  % "folder" - folder name
  % for actual loading uses loadimageseries routine
  global rawPhaseData rawPhaseFolder rawS1Data rawS2Data imsizes imageFolders imageLimits %#ok<NUSED>
  if n==1, str='rawPhaseData';rawPhaseFolder = folder; end
  % if n==2, str='rawFMData'; end
  if n==3, str='rawS1Data'; end
  if n==4, str='rawS2Data'; end
  filenames = ''; %#ok<NASGU>
  dirname = ''; %#ok<NASGU>
  imageFolders{n} = folder;
  eval(['[' str ', filenames, folder] = loadimageseries(folder,1);']);
  eval(['cls = class(' str ');']);
  if strcmp(cls,'uint8'), lng=8; elseif strcmp(cls,'uint16'), lng=16; elseif strcmp(cls,'uint32'), lng=32;
    else gdisp('No images loaded');return; end
  eval(['imageLimits{n} = 2^' num2str(lng) '*mean(stretchlim(' str ',[0.0001 0.9999]),2);']);
  %eval(['imageLimits{n} = im2double([min(min(min(' str '))) max(max(max(' str ')))]);']);
  if(folder==-1), return; end; 
  eval(['imsizes(n,:) = [size(' str ',1) size(' str ',3) size(' str ',3)];']);
  gdisp(['Loaded ' num2str(imsizes(n,3)) ' files'])
  updateimsizes
end

function updateimsizes
  % Updates the structure "imsizes" chat contains the information about 
  % the size of each images stack (phase, extra, signa1, signal2) and
  % the area occupied by the meshes. Needed for display purposes.
  global rawPhaseData rawS1Data rawS2Data imsizes cellList regionSelectionRect
  imsizesold = imsizes;
  imsizes(1,:)=[size(rawPhaseData,1) size(rawPhaseData,2) size(rawPhaseData,3)];
  % imsizes(2,:)=[size(rawFMData,1) size(rawFMData,2) size(rawFMData,3)];
  imsizes(3,:)=[size(rawS1Data,1) size(rawS1Data,2) size(rawS1Data,3)];
  imsizes(4,:)=[size(rawS2Data,1) size(rawS2Data,2) size(rawS2Data,3)];
  if ~isempty(cellList)
    xmax=0;
    ymax=0;
    frmmax = length(cellList);
    for frm=1:max(1,floor(frmmax/3)):frmmax
      for i=1:length(cellList{frm})
        if ~isempty(cellList{frm}{i})
          box = cellList{frm}{i}.box;
          xmax = max(xmax,box(2)+box(4));
          ymax = max(ymax,box(1)+box(3));
        end
      end
    end
    imsizes(end,:) = max([imsizes(1:end-1,:);[xmax ymax length(cellList)]]);
  else
    imsizes(end,:) = max(imsizes(1:end-1,:));
  end
  if imsizes(end,1)==0, imsizes(end,:) = [400 500 1]; end
  if ~prod((imsizesold(end,:)==imsizes(end,:))+0)
    regionSelectionRect = [];
  end
end

%% Parameters global functions

function res = getparamstring(hndls)
  str = get(hndls.params,'String');
  if ~iscell(str) && size(str,1)==1
    res = textscan(str,'%s','delimiter','');
  elseif ~iscell(str) && size(str,1)>1
    for i=1:size(str,1)
      tmp1 = strtrim(str(i,:));
      if ~isempty(tmp1)
        tmp2 = textscan(tmp1,'%s','delimiter','');
        str2{i,1} = tmp2{1}{1};
      end
    end
    res = str2;
  else
    res = str;
  end
  if length(res)==1
    res=res{1};
  end
end

function res=loadparam(filename)
  % loads parameters either from parameters file or from meshes file, in
  % first case the file is a text file, in the second - a string variable
  % the function sends the data to "parseparams" function (updates the 
  % global variable "p") and returns the data string to display
  if length(filename)>4 && strcmp(filename(end-3:end),'.set')
    try
      res = fileread(filename); % Saving/loading format
    catch
      errordlg(['Could not open parameters file! Make sure the file exists and'...
        ' the parameters are saved in the correct format.']);
      return;
    end;
  elseif length(filename)>4 && strcmp(filename(end-3:end),'.mat')
    res = load(filename);
    if isempty(res)
      errordlg('Could not parce parameters! No parameters saved in this file.');
      return;
    else
      try
        res = {res.paramString};
        % res = textscan(str,'%s','delimiter','');
      catch
        errordlg('Could not parse parameters! Make sure the parameters are saved in the correct format.');
        return;
      end
    end
  else
    errordlg('Could not open parameters file! File extension must be ".set" or ".mat".');
    return;
  end
  if ~iscell(res)
    res=textscan2(res);
    % res = textscan(res,'%s','delimiter','');
    % res=res{1};
  end
  parseparams(res)
end

function [pOut] = returnparam()
global p;
pOut = p;
end

function res=textscan2(str)
  str = strtrim(str);
  str = regexprep(str,char([13 10]),char(10));
  str = regexprep(str,char([92 110]),char(10));
  pos = [0 sort([strfind(str,char(10)) strfind(str,'\n')]) length(str)+1];
  res = {};
  for i=1:length(pos)-1
    if pos(i)+1<=pos(i+1)-1
      res{i,1} = str(pos(i)+1:pos(i+1)-1);
    else
      res{i,1} = '';
    end
  end
end

function parseparams(str)
  % takes parameters string (loaded or modified by the user) and updates
  % the global variable "p" - parameters structure
  global p handles
  names = {};
  values = {};
  if length(str)==1, str = str{1}; end
  for i=1:length(str)
    res = textscan([char(str{i}) ' '], '%s%s', 'commentstyle', '%','delimiter','=');
    if ~isempty(res{1})
      names = [names res{1}];
      values = [values res{2}];
    end
  end
  try
    p=[];
    for i=1:length(names)
      if ~isempty(str2num( values{i}))
        eval(['p.' names{i} '=[' values{i} '];']);
      else
        eval(['p.' names{i} '=[' char(39) strtrim(values{i}) char(39) '];']);
      end
    end
    % Temporarily added for new parameters
    if ~isfield(p,'aligndepth'), p.aligndepth = 1; end
    if ~isfield(p,'getmesh'), p.getmesh = true; end
    if ~isfield(p,'meshStep'), p.meshStep = 1; end
    if ~isfield(p,'fmeshstep'), p.fmeshstep = p.meshStep; end
    if ~isfield(p,'scaleFactor'), p.scaleFactor = 1; end
    if ~isfield(p,'meshTolerance'), p.meshTolerance = 0.01; end
    if ~isfield(p,'meshWidth') && ~isfield(p,'cellwidth'), p.meshWidth = 12; end
    if ~isfield(p,'meshWidth') && isfield(p,'cellwidth'), p.meshWidth = p.cellwidth*1.5; end % !!!
    if ~isfield(p,'neighRepA'), p.neighRepA = 0; end
    if ~isfield(p,'useExtraData'), p.useExtraData = false; end
    if ~isfield(p,'invertimage'), p.invertimage = p.useExtraData; end
    if ~isfield(p,'maxRegNumber'), p.maxRegNumber = 100000; end
    if ~isfield(p,'joindilate'), p.joindilate = 1; end
    if ~isfield(p,'edgemode'), if isfield(p,'logmode'), p.edgemode = p.logmode; else p.edgemode = 'log'; end; end
    if isfield(p,'edgedetection') && ~p.edgedetection, p.edgemode = 'none'; end
    if ~isfield(p,'edgeSigmaL'), if isfield(p,'edgeSigma'), p.edgeSigmaL = p.edgeSigma; else p.edgeSigmaL = 1; end; end
    if ~isfield(p,'edgeSigmaV'), p.edgeSigmaV = 0.5; end
    if ~isfield(p,'valleythresh1'), if isfield(p,'valleythres1'), p.valleythresh1 = p.valleythres1*300; else p.valleythresh1 = 0; end; end
    if ~isfield(p,'valleythresh2'), if isfield(p,'valleythres2'), p.valleythresh2 = p.valleythres2*300; else p.valleythresh2 = 1; end; end
    if ~isfield(p,'fitDisplay'), p.fitDisplay = false; end
    if ~isfield(p,'fitDisplay1'), p.fitDisplay1 = false; end
    if ~isfield(p,'opennum'), p.opennum = 0; end
    if ~isfield(p,'threshminlevel'), p.threshminlevel = 0; end
    if ~isfield(p,'approxSignal'), p.approxSignal = 0; end
    if ~isfield(p,'forceindframes'), p.forceindframes = 0; end
    if ~isfield(p,'preProcSzObj'),p.preProcSzObj = 150; end
    if ~isfield(p,'preProcWienerBox'),p.preProcWienerBox = [5,5] ; end
    catch ME
      error('The format of one or more parameters is incorrect or parameters missing');
    return;
  end
end

%% Saving mesh global functions

function param = loadmesh(filename)
  % loads mesh from the specified file, updates "cellList" and
  % "cellListN" structures, sets "selectedList" to empty = no cells
  % selected, returns parameters structure (which can be used or not,
  % depending on the "Load params" checkbox)
  global cellList cellListN selectedList
  warning off
  l=load(filename,'cellList','cellListN','paramString');
  warning on
  cellList = l.cellList;
  cellListN = l.cellListN;
  selectedList = [];
  if isfield(l,'paramString')
    param = l.paramString;
  else
    param = [];
  end
  gdisp(['Meshes loaded from file ' filename])
end

function savemesh(filename,lst,savelst,range)
  % saves meshes:
  % <lst> - list of cells (only used if <savelst>=true)
  % <range> - 2-element array indicating the frame range to save
  global p coefPCA weights mCell rawPhaseFolder cellList cellListN handles shiftframes shiftfluo %#ok<NUSED>
  if isempty(filename) || isempty(cellList), return; end
  if isempty(range) || ~isnumeric(range), range = [1 length(cellList)]; end
  if length(range)==1, range=[range range]; end
  range = [max(range(1),1),min(range(2),length(cellList))];
  if range(2)-range(1)<0, return; end
  %paramString = getparamstring(handles);
  if (isempty(lst) || ~savelst) && range(1)==1 && range(2)==length(cellList) % Saving whole cellList
    if ~isempty(dir(filename)), delete(filename); pause(0.1); end
    %save(filename,'cellList','cellListN','p','coefPCA','weights','mCell','rawPhaseFolder','paramString','shiftframes','shiftfluo');
    save(filename,'cellList','cellListN','p','coefPCA','weights','mCell','rawPhaseFolder','shiftframes','shiftfluo');
  elseif (isempty(lst) || ~savelst) % Saving all cells in a range of frames
    cellListTmp = {};
    for j=range(1):range(2)
      cellListTmp{j} = cellList{j};
    end
    if ~isempty(dir(filename)), delete(filename); pause(0.1); end
    %save(filename,'cellListN','p','coefPCA','weights','mCell','rawPhaseFolder','paramString','shiftframes','shiftfluo');
    save(filename,'cellListN','p','coefPCA','weights','mCell','rawPhaseFolder','shiftframes','shiftfluo');
    savetmp(filename,cellListTmp)
  else % Saving selected cells in a range of frames
    cellListTmp = {};
    for j=range(1):range(2)
      if j~=range(1), lst = selNewFrame(lst,j-1,j); end
      for i=lst
        cellListTmp{j}{i} = cellList{j}{i};
      end
    end
    if ~isempty(dir(filename)), delete(filename); pause(0.1); end
    %save(filename,'cellListN','p','coefPCA','weights','mCell','rawPhaseFolder','paramString','shiftframes','shiftfluo');
    save(filename,'cellListN','p','coefPCA','weights','mCell','rawPhaseFolder','shiftframes','shiftfluo');
    savetmp(filename,cellListTmp)
  end
  % display if succeded
  [tmp,filename,ext] = fileparts(filename);
  if isempty(ext), ext='.mat'; end
  gdisp(['Analysis saved to file ' filename ext])
end

function savetmp(filename,cellList)
  % This function adds cellList to a file
  % to use with 'savemesh' function only
  save(filename,'-append','cellList')
end

%% Computing ancestors/descendants global functions

function newlst2 = selNewFrame(lst,oldframe,newframe)
  % updates the list of selected cells ("lst" as selected on frame 
  % "oldframe") and oputputs as mewlst2 (on frame "newframe")
  % (updated in version 0.933)
  global cellList cellListN
  if oldframe==newframe, newlst2=lst; return; end
  newlst2 = [];
  if max(newframe,oldframe)>length(cellList) || max(newframe,oldframe)>length(cellListN) ||...
      isempty(lst) || isempty(cellList{oldframe}) || isempty(cellList{newframe}) ||...
      cellListN(oldframe)~=cellListN(newframe), 
    return
  end
  for cell=1:length(cellList{oldframe})
    if ~isempty(cellList{oldframe}{cell}) && isfield(cellList{oldframe}{cell},'timelapse')
      if cellList{oldframe}{cell}.timelapse
        break % is the first cell has timelapse=1 property, proceed with selection
      else
        return % is the first cell has timelapse=0 property, no cells should be selected
      end
    end
  end
  if newframe>oldframe
    maxcell = cellListN(oldframe);
    for i=lst
      if ~isempty(cellList{oldframe}{i}) && isfield(cellList{oldframe}{i},'ancestors')
        if i<=length(cellList{newframe}) && ~isempty(cellList{newframe}{i}) && ...
            isfield(cellList{oldframe}{i},'birthframe') && isfield(cellList{newframe}{i},'birthframe') && ...
            cellList{newframe}{i}.birthframe ~= cellList{oldframe}{i}.birthframe
          newlst2 = []; % birth frames does not match => not timelapse or discontinued lineage
          return
        end
        igen = length(cellList{oldframe}{i}.ancestors)+length(cellList{oldframe}{i}.descendants);
        maxcelli = maxcell*2^igen;
        newlst3 = [];
        for j=i:maxcelli:length(cellList{newframe})
          if (~isempty(cellList{newframe}{j}) && (i==j || ismember(i,cellList{newframe}{j}.ancestors))) ...
              || ((newframe==(oldframe+1)) && (j==i)) % this line is for processing selected cells
            newlst3 = [newlst3 j];
          end
        end
        newlst2 = [newlst2 newlst3];
      end
    end
  else
    newlst = lst;
    for i=newlst % keep only the cells existing on the new frame
      if newframe<=length(cellList) && i<=length(cellList{newframe}) && ~isempty(cellList{newframe}{i})
        if isfield(cellList{oldframe}{i},'birthframe') && isfield(cellList{newframe}{i},'birthframe') && ...
            cellList{newframe}{i}.birthframe ~= cellList{oldframe}{i}.birthframe
          newlst2 = []; % birth frames does not match => not timelapse or discontinued lineage
          return
        end
        newlst2 = [newlst2 i];
      end
    end
  end
end

function res = getdaughter(cell,generation,MaxCell)
  % gets the list of a given cell progeny born after the starting plane
  % in the indicated generation and for the indicated starting maximum
  % number of cells ("MaxCell", number of cells on the first frame)
  res = cell + MaxCell * 2^(max(1,generation)-1+ceil(log2(ceil(cell/MaxCell))));
end

function res = getmother(cell,MaxCell)
  if cell>MaxCell
    MaxCell2 = MaxCell * 2^(-1+ceil(log2(ceil(cell/MaxCell))));
    res = cell - MaxCell2;
  else
    res = [];
  end
end

function res = getallchildren(cell,MaxCell,cellsOnFrame)
  % gets all POSSIBLE progeny of "cell" on the current frame
  % (the real progeny is computed in the parent function "selNewFrame")
  MaxCell2 = MaxCell * 2^(ceil(log2(ceil(cell/MaxCell))));
  res = cell:MaxCell2:cellsOnFrame;
end

% --- End of computing ancestors/descendants global functions ---

%% Adding signal global functions

function addSignalToFrame(frame,inchannel,outchannel,proccells,rsz,apr,shiftfluo)
  global cellList rawPhaseData rawS1Data rawS2Data
  if length(cellList)<frame || isequal(outchannel,'-'), return; end
  if frame<1, frame=1; end
  if exist('aip','file')==0, apr=true; end
  sf = [0 0];
  switch inchannel
   case 1
     if ~isempty(who('rawPhaseData'))
       if ~isempty(rawPhaseData)
         if frame<=size(rawPhaseData,3)
          img = im2double(rawPhaseData(:,:,frame));
          img = max(max(img))-img;
         end
       end
     end
  %  case 2
  %    if ~isempty(who('rawFMData'))
  %      if ~isempty(rawFMData)
  %        if frame<=size(rawFMData,3)
  %          img = im2double(rawFMData(:,:,frame));
  %        end
  %      end
  %    end
   case 3
     if ~isempty(who('rawS1Data'))
       if ~isempty(rawS1Data)
         if frame<=size(rawS1Data,3)
          if ~isempty(shiftfluo), sf = shiftfluo(1,:); end
          img = im2double(rawS1Data(:,:,frame));
         end
       end
     end
   case 4
     if ~isempty(who('rawS2Data'))
       if ~isempty(rawS2Data)
         if frame<=size(rawS2Data,3)
          if ~isempty(shiftfluo), sf = shiftfluo(2,:); end
          img = im2double(rawS2Data(:,:,frame));
         end
       end
     end
   otherwise
    return
  end
  if isempty(who('img')), return; end
  if isempty(proccells), proccells = 1:length(cellList{frame}); end
  proccells2 = [];
  for i=proccells
    if length(cellList{frame})>=i
      if ~isempty(cellList{frame}{i})
        proccells2 = [proccells2 i];
      end
    end
  end
  proccells = proccells2;
  for cell = proccells
    if isfield(cellList{frame}{cell},'mesh') && size(cellList{frame}{cell}.mesh,1)>1 && apr
      mesh = cellList{frame}{cell}.mesh;
      mesh(:,[1 3])=mesh(:,[1 3])+sf(1); mesh(:,[2 4])=mesh(:,[2 4])+sf(2);
      S = getOneSignalM(mesh,cellList{frame}{cell}.box,img,rsz);
    elseif isfield(cellList{frame}{cell},'mesh') && size(cellList{frame}{cell}.mesh,1)>1 && ~apr
      mesh = cellList{frame}{cell}.mesh;
      mesh(:,[1 3])=mesh(:,[1 3])+sf(1); mesh(:,[2 4])=mesh(:,[2 4])+sf(2);
      S = getOneSignalC(mesh,cellList{frame}{cell}.box,img,rsz);
    elseif isfield(cellList{frame}{cell},'contour') && apr
      contour = cellList{frame}{cell}.contour;
      contour(:,1)=contour(:,1)+sf(1); contour(:,2)=contour(:,2)+sf(2);
      S = getOneSignalContourM(contour,cellList{frame}{cell}.box,img,rsz);
    elseif isfield(cellList{frame}{cell},'contour') && ~apr
      contour = cellList{frame}{cell}.contour;
      contour(:,1)=contour(:,1)+sf(1); contour(:,2)=contour(:,2)+sf(2);
      S = getOneSignalContourC(contour,cellList{frame}{cell}.box,img,rsz);
    else
      S = [];
    end
    if isnumeric(outchannel)
      eval(['cellList{frame}{cell}.signal' num2str(outchannel) '=S;']);
    else
      eval(['cellList{frame}{cell}.' outchannel '=S;']);
    end
  end
end

function sgn = getOneSignalM(mesh,mbox,img,rsz)
  % This function integrates signal confined within each segment of the
  % mesh. Two versions are provided. This one is a faster but less
  % precise MATLAB-based approximation.
  % 
  sgn = [];
  if length(mesh)>1
    img2 = imcrop(img,mbox);
    if rsz>1
      img2 = imresize(img2,rsz);
    end
    for i=1:size(mesh,1)-1
      plgx = rsz*([mesh(i,[1 3]) mesh(i+1,[3 1])]-mbox(1)+1);
      plgy = rsz*([mesh(i,[2 4]) mesh(i+1,[4 2])]-mbox(2)+1);
      mask = poly2mask(plgx,plgy,size(img2,1),size(img2,2));
      s = sum(sum(mask));
      if s>0, s=(rsz^(-2))*sum(sum(mask.*img2))*polyarea(plgx,plgy)/s; end
      sgn = [sgn;s];
    end
  end
end

function sgn = getOneSignalC(mesh,mbox,img,rsz)
  % This function integrates signal confined within each segment of the
  % mesh. Two versions are provided. This one is a slower but more
  % precise C-based function (requires 'aip' function to run).
  %
  sgn = zeros(size(mesh,1)-1,1);
  if length(mesh)>1
    img2 = imcrop(img,mbox);
    if rsz>1
      img2 = imresize(img2,rsz);
    end
    c = repmat(1:size(img2,2),size(img2,1),1);
    r = repmat((1:size(img2,1))',1,size(img2,2));
    icw = isContourClockwise([mesh(:,1);flipud(mesh(2:end-1,3))],[mesh(:,2);flipud(mesh(2:end-1,4))]);
    maxi = size(mesh,1)-1;
    for i=1:maxi
      plgx = rsz*([mesh(i,[1 3]) mesh(i+1,[3 1])]-mbox(1)+1)';
      plgy = rsz*([mesh(i,[2 4]) mesh(i+1,[4 2])]-mbox(2)+1)';
      if icw, plgx = flipud(plgx); plgy = flipud(plgy); end % making contours clockwise
      [plgx2,plgy2] = expandpoly(plgx,plgy,1.42,1);
      if sum(isnan(plgx2))>0, sgn=[]; return; end
      mask = poly2mask(plgx2,plgy2,size(img2,1),size(img2,2));
      f = find(mask);
      int = 0;
      if ~isempty(f)
        for px=f'
          if (i==1||i==maxi) && ~isContourClockwise(plgx,plgy)
            plgx = flipud(plgx); plgy = flipud(plgy);
            int = int-img2(px)*aip(plgx,plgy,c(px)-0.5,r(px)-0.5);
          else
            int = int+img2(px)*aip(plgx,plgy,c(px)-0.5,r(px)-0.5);
          end
        end
      end
      sgn(i) = sum(int);
    end
  end
end

function sgn = getOneSignalContourM(contour,mbox,img,rsz)
  % This function integrates signal confined within the contour (which is 
  % not segmented in the way the mesh is). This is a faster but less
  % precise MATLAB-based version.
  % 
  if length(contour)>1
    img2 = imcrop(img,mbox);
    if rsz>1
      img2 = imresize(img2,rsz);
    end
    plgx = rsz*(contour(:,1)-mbox(1)+1);
    plgy = rsz*(contour(:,2)-mbox(2)+1);
    mask = poly2mask(plgx,plgy,size(img2,1),size(img2,2));
    sgn = sum(sum(mask));
    if sgn>0, sgn=(rsz^(-2))*sum(sum(mask.*img2))*polyarea(plgx,plgy)/sgn; end
  end
end

function sgn = getOneSignalContourC(contour,mbox,img,rsz)
  % This function integrates signal confined within the contour (which is 
  % not segmented in the way the mesh is). This is a slower but more
  % precise C-based function (requires 'aip' function to run).
  % 
  if length(contour)>1
    img2 = imcrop(img,mbox);
    if rsz>1
      img2 = imresize(img2,rsz);
    end
    c = repmat(1:size(img2,2),size(img2,1),1);
    r = repmat((1:size(img2,1))',1,size(img2,2));
    plgx = rsz*(contour(:,1)-mbox(1)+1);
    plgy = rsz*(contour(:,2)-mbox(2)+1);
    [plgx,plgy] = poly2cw(plgx,plgy); % making contour clockwise
    [plgx2,plgy2] = expandpoly(plgx,plgy,1.42,1); % making contours clockwise
    if sum(isnan(plgx2))>0, sgn=[]; return; end
    mask = poly2mask(plgx2,plgy2,size(img2,1),size(img2,2));
    f = find(mask);
    int = 0;
    if ~isempty(f)
      for px=f'
         int = int+img2(px)*aip(plgx,plgy,c(px)-0.5,r(px)-0.5);
      end
    end
    sgn = sum(int);
  end
end
  
% --- End of adding signal global functions ---

% -------------------------------------------------------------------------
%% ---- Processing functions ----

function proccells = processFrame(frame,proccells)
% This function is for processing all frames except the first one
% 
% It is only looking for the cells that existed before and splits them if
% necessary.
%SH150225 modified to have the same changes as processFrameI

global maskdx maskdy imsizes cellList cellListN rawPhaseData se p

cellListN(frame) = cellListN(frame-1);
if p.invertimage
  if frame>size(rawPhaseData,3), return; end
  img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
else
  if frame>size(rawPhaseData,3), return; end
  img = rawPhaseData(:,:,frame);
end
img = preProcessPhIm(img,p.preProcSzObj,p.preProcWienerBox);
imge = img2imge(img,p.erodeNum);
imge16 = img2imge16(img,p.erodeNum);
thres = graythreshreg(imge,p.threshminlevel);
bgr = phasebgr(imge,thres);

allMap = zeros(size(img));
for cell = 1:length(cellList{frame-1})
  if isempty(cellList{frame-1}{cell}), continue; end
  mesh = cellList{frame-1}{cell}.mesh;
  if size(mesh,2)==4
    allMap = allMap +... 
      imerode(roipoly(allMap,[mesh(:,1);flipud(mesh(:,3))],[mesh(:,2);flipud(mesh(:,4))]),se);
  end
end
allMap = min(allMap,1);

[extDx,extDy] = getExtForces(imge,imge16,p);
if isempty(extDx), gdisp('Processing timelapse frame failed: unable to get energy'); return; end

%cell = 1;
%cellmax = length(cellList{frame-1});
if isempty(proccells)
  for i=1:length(cellList{frame-1})
    if ~isempty(cellList{frame-1}{i})
      proccells = [proccells i];
    end
  end
end
%cellmaxCell = cellmax+1;

% Check 2
for cell = proccells % parfor
  if cell>length(cellList{frame-1}) || isempty(cellList{frame-1}{cell}), continue; end
  %gdisp(['processing cell ' num2str(cell)])
  
  % get previous frame data
  prevStruct = cellList{frame-1}{cell};
  if size(prevStruct.model,2)==1 || size(prevStruct.model,2)==1 
    pcCell = reshape(prevStruct.model,[],1);
  else
    pcCell = prevStruct.model;
  end
  cCell = model2geom(pcCell,p.algorithm);
  
  if ~isempty(cCell)
    roiBox(1:2) = round(max(min(cCell(:,1:2))-p.roiBorder,1));
    roiBox(3:4) = min(round(max(cCell(:,1:2))+p.roiBorder),[size(img,2) size(img,1)])-roiBox(1:2);
    if min(roiBox(3:4))<=0, cCell=[]; end
  end
  if ~isempty(cCell)
    pcCellRoi = model2box(pcCell,roiBox,p.algorithm);

    % crop the image and the energy/force maps
    roiImg = imcrop(imge,roiBox);
    roiExtDx = imcrop(extDx,roiBox);
    roiExtDy = imcrop(extDy,roiBox);
    roiAmap = imcrop(allMap,roiBox);

    % build pmap (for cell repulsion)
    mesh = cellList{frame-1}{cell}.mesh;
    if size(mesh,2)==4
      roiAmap = max(0,roiAmap - ... 
        imdilate(roipoly(roiAmap,[mesh(:,1);flipud(mesh(:,3))]-roiBox(1),...
                     [mesh(:,2);flipud(mesh(:,4))]-roiBox(2)),se));
    end
    pmap = roiAmap;
    f1=true;
    while f1
      pmap1 = roiAmap + imerode(pmap,se);
      f1 = max(max(pmap1-pmap))>0;
      pmap = pmap1;
    end;
    pmapDx = imfilter(pmap,maskdx,'replicate'); % distance forces
    pmapDy = imfilter(pmap,maskdy,'replicate'); 
    roiExtDx = roiExtDx + p.neighRep*pmapDx;
    roiExtDy = roiExtDy + p.neighRep*pmapDy;

    if ismember(p.algorithm,[2 3])
      [pcCellRoi,fitquality] = align(roiImg,roiExtDx,roiExtDy,roiAmap,pcCellRoi,p,false,roiBox,thres,[frame cell]);
    elseif p.algorithm == 4
      [pcCellRoi,fitquality] = align4(roiImg,roiExtDx,roiExtDy,roiAmap,pcCellRoi,p,roiBox,thres,[frame cell]);
    end
    % TEST HERE
    %gdisp(['fitquality aligning = ' num2str(fitquality)])

    % obtaining the shape of the cell in geometrical representation
    pcCellBox = pcCell;%SH150225
    pcCell = box2model(pcCellRoi,roiBox,p.algorithm);
    cCell = model2geom(pcCell,p.algorithm);

    %try splitting
    if p.algorithm==4, isct=0; else isct = intersections(cCell); end
    cellarea = 0;
    if ~isempty(isct) && ~isct
      mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth);
      if length(mesh)>1
        cellarea = polyarea([mesh(:,1);flipud(mesh(:,3))],[mesh(:,2);flipud(mesh(:,4))]);
        roiMesh = mesh - repmat(roiBox([1 2 1 2])-1,size(mesh,1),1);
        res=isDivided(roiMesh,roiImg,p.splitThreshold,bgr);
        if res==-1, mesh=0; end
      end
    end
  end
  
  % delete all descendants of the cell if already present
  if ~isfield(p,'delchildrenonredetect') || ~p.delchildrenonredetect
    fmax = length(cellList);
  else
    fmax = frame;
  end
  for frm = frame:fmax
    dlst = selNewFrame(cell,frame-1,frm);
    for cl = dlst;
      if cl<=length(cellList{frm}) && ~isempty(cellList{frm}{cl})
        cellList{frm}{cl} = [];
      end
    end
  end

  % checking quality and storing the cell
  if ~isempty(cCell) && ~isempty(isct) && ~isct && fitquality<p.fitqualitymax && ...
       min(cCell(:,1))>= p.imRoi(1) &&...
       min(cCell(:,2))>= p.imRoi(2) &&...
       max(cCell(:,1))<=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)) && ...
       max(cCell(:,2))<=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)) && ...
       length(mesh)>1 && p.areaMin<cellarea && p.areaMax>cellarea

      %calculate the error estimates SH150225
      %NOT IMPLEMENTED YET - because roiMask not calculated
      %errorEstimate = getMTrackErrorEstimate(p,pcCellBox,mesh,fitquality,roiImg,roiExtDx,roiExtDy,roiMask,roiBox);
      
    % if the cell passed the quality test and it is not on the
    % boundary of the image - store it
    if res>0 % Splitting the cell
      % mesh1 = flipud(mesh(1:res-1,:)); % daughter cell
      % mesh2 = mesh(res+1:end,:);
      % select the cells so that the mesh always goes from the stalk pole
      mesh1 = flipud(mesh(res+1:end,:)); % daughter cell
      mesh2 = mesh(1:res-1,:); % mother cell

      if ismember(p.algorithm,[2 3])
        pcCell1 = splitted2model(mesh1,p);
        pcCell1 = model2box(pcCell1,roiBox,p.algorithm);
        pcCell1 = align(roiImg,roiExtDx,roiExtDy,roiAmap,pcCell1,p,false,roiBox,thres,[frame cell]);
        pcCell1 = box2model(pcCell1,roiBox,p.algorithm);
        cCell1 = model2geom(pcCell1,p.algorithm);
        pcCell2 = splitted2model(mesh2,p);
        pcCell2 = model2box(pcCell2,roiBox,p.algorithm);
        pcCell2 = align(roiImg,roiExtDx,roiExtDy,roiAmap,pcCell2,p,false,roiBox,thres,[frame cell]);
        pcCell2 = box2model(pcCell2,roiBox,p.algorithm);
        cCell2 = model2geom(pcCell2,p.algorithm);
      else
        pcCell1 = align4IM(mesh1,p);
        pcCell1 = model2box(pcCell1,roiBox,p.algorithm);
        pcCell1 = align4(roiImg,roiExtDx,roiExtDy,roiAmap,pcCell1,p,roiBox,thres,[frame cell]);
        pcCell1 = box2model(pcCell1,roiBox,p.algorithm);
        cCell1 = pcCell1;
        pcCell2 = align4IM(mesh2,p);
        pcCell2 = model2box(pcCell2,roiBox,p.algorithm);
        pcCell2 = align4(roiImg,roiExtDx,roiExtDy,roiAmap,pcCell2,p,roiBox,thres,[frame cell]);
        pcCell2 = box2model(pcCell2,roiBox,p.algorithm);
        cCell2 = pcCell2;
      end
      
      mesh1 = model2mesh(cCell1,p.meshStep,p.meshTolerance,p.meshWidth);
      mesh2 = model2mesh(cCell2,p.meshStep,p.meshTolerance,p.meshWidth);
      cellStruct1.algorithm = p.algorithm;
      cellStruct2.algorithm = p.algorithm;
      cellStruct1.birthframe = frame;
      cellStruct2.birthframe = prevStruct.birthframe;
      if size(pcCell1,2)==1, model=pcCell1'; else model=pcCell1; end
      cellStruct1.model = model;
      if size(pcCell2,2)==1, model=pcCell2'; else model=pcCell2; end
      cellStruct2.model = model;
      cellStruct1.polarity = 1;
      cellStruct2.polarity = 1;
      cellStruct1.mesh = mesh1;
      cellStruct2.mesh = mesh2;
      cellStruct1.stage = 1;
      cellStruct2.stage = 1;
      cellStruct1.timelapse = 1;
      cellStruct2.timelapse = 1;
      cellStruct1.divisions = []; % frame?
      cellStruct2.divisions = [prevStruct.divisions frame];
      cellStruct1.box = roiBox;
      cellStruct2.box = roiBox;
      daughter = getdaughter(cell,length(cellStruct2.divisions),cellListN(frame));
      cellStruct1.ancestors = [prevStruct.ancestors cell];
      cellStruct2.ancestors = prevStruct.ancestors;
      cellStruct1.descendants = [];
      cellStruct2.descendants = [prevStruct.descendants daughter];
      %cellStruct1.errorEstimate = errorEstimate;
      %cellStruct2.errorEstimate = errorEstimate;
      cellList{frame}{cell} = cellStruct2; % mother cell keeps the number
      if (frame>1)&&(length(cellList{frame-1})>=cell)&&~isempty(cellList{frame-1}{cell}), cellList{frame-1}{cell}.timelapse=1; end
      if daughter<p.maxCellNumber
        cellList{frame}{daughter} = cellStruct1;
        proccells = [proccells daughter];
      else
        gdisp(['cell ' num2str(daughter) ' born still because of overpopulation!'])
      end
      gdisp(['cell ' num2str(cell) ' splitted, cell ' num2str(daughter) ' was born'])
      continue;
    end
    
    cellStruct.birthframe = prevStruct.birthframe;
    cellStruct.algorithm = p.algorithm;
    if size(pcCell,2)==1, model=pcCell'; else model=pcCell; end
    cellStruct.model = model;
    cellStruct.mesh = mesh;
    cellStruct.polarity = prevStruct.polarity;
    cellStruct.ancestors = prevStruct.ancestors;
    cellStruct.descendants = prevStruct.descendants;
    cellStruct.stage = 1;
    cellStruct.timelapse = 1;
    %cellStruct.stage = getStage(roiImg,cellStruc.model);
    cellStruct.divisions = prevStruct.divisions;
    cellStruct.box = roiBox;
    %cellStruct.errorEstimate = errorEstimate;
    cellList{frame}{cell} = cellStruct;
    if (frame>1)&&(length(cellList{frame-1})>=cell)&&~isempty(cellList{frame-1}{cell}), cellList{frame-1}{cell}.timelapse=1; end
    gdisp(['fitting cell ' num2str(cell) ' - passed and saved'])
    continue;
  else
    cellStruct = prevStruct;
    if size(pcCell,2)==1, model=pcCell'; else model=pcCell; end
    cellStruct.model = model;
    cellStruct.mesh = 0;
    cellStruct.stage = 0;
    cellList{frame}{cell} = cellStruct; % remove this line to not save empty meshes
    % if the cell is not on the image border OR it did not pass the
    % quality test - split it
    reason = 'unknown';
    if isempty(pcCell) || isempty(cCell), reason = 'no cell found';
    elseif fitquality>=p.fitqualitymax, reason = 'bad fit quality';
    elseif min(cCell(:,1))<=p.imRoi(1), reason = 'cell on x=0 boundary';
    elseif min(cCell(:,2))<=p.imRoi(2), reason = 'cell on y=0 boundary'; 
    elseif max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)), reason = 'cell on x=max boundary'; 
    elseif max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)), reason = 'cell on y=max boundary'; 
    elseif isct, reason = 'model has intersections';
    elseif length(mesh)<=1 || cellarea==0, reason = 'problem getting mesh';
    elseif p.areaMin>=cellarea, reason = 'cell too small';
    elseif p.areaMax<=cellarea, reason = 'cell too big';
    end
    gdisp(['fitting cell ' num2str(cell) ' - quality check failed - ' reason]);  
    
  end
end
end % processFrame function



function testModeSegment(frame,processregion)
% Display of the segmented image for parameter testing purposes

  global rawPhaseData p;

  if checkparam(p,'invertimage','algorithm','erodeNum','thresFactorM','opennum','edgeSigmaL','edgeSigmaV','valleythresh1','valleythresh2')
    gdisp('Segmentation failed: one or more required parameters not provided.');
    return
  end 
  if frame>size(rawPhaseData,3), gdisp('Segmentation failed: no phase images loaded'); return; end
  if p.invertimage
    img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
  else
    img = rawPhaseData(:,:,frame);
  end
   img = preProcessPhIm(img,p.preProcSzObj,p.preProcWienerBox);
  imge = img2imge(img,p.erodeNum);
  imge16 = img2imge16(img,p.erodeNum);
  
  thres = graythreshreg(imge,p.threshminlevel);
  regions0 = getRegions(p.edgemode,imge,thres,p.thresFactorM,p.opennum,imge16,p.edgeSigmaL,p.edgeSigmaV,p.valleythresh1,p.valleythresh2,p.imRoi);
  if ~isempty(processregion)
    crp = regions0(processregion(2)+(0:processregion(4)),processregion(1)+(0:processregion(3)));
    regions0 = regions0*0;
    regions0(processregion(2)+(0:processregion(4)),processregion(1)+(0:processregion(3))) = crp;
    regions0 = bwlabel(regions0>0,4);
  end
  
  figure
  if isempty(processregion)
    imshow(regions0>0,[]);set(gca,'pos',[0 0 1 1])
  else
    imshow(imcrop(regions0,processregion)>0,[]);set(gca,'pos',[0 0 1 1])
  end
end


function processFrameI(frame,tl,processregion)
% This function is for processing the first frame only
% 
% It searches for the cells on the frame and then does refinement

global imsizes se maskdx maskdy cellList cellListN rawPhaseData p

if p.invertimage
  if frame>size(rawPhaseData,3), return; end
  img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
else
  if frame>size(rawPhaseData,3), return; end
  img = rawPhaseData(:,:,frame);
end
img = preProcessPhIm(img,p.preProcSzObj,p.preProcWienerBox);
imge = img2imge(img,p.erodeNum);
imge16 = img2imge16(img,p.erodeNum);
imgemax = max(max(imge));
thres = graythreshreg(imge,p.threshminlevel);
bgr = phasebgr(imge,thres);

regions0 = getRegions(p.edgemode,imge,thres,p.thresFactorM,p.opennum,imge16,p.edgeSigmaL,p.edgeSigmaV,p.valleythresh1,p.valleythresh2,p.imRoi);

[extDx,extDy] = getExtForces(imge,imge16,p);
if isempty(extDx), gdisp('Processing independent frame failed: unable to get energy'); return; end

if ~isempty(processregion)
  crp = regions0(processregion(2)+(0:processregion(4)),processregion(1)+(0:processregion(3)));
  regions0 = regions0*0;
  regions0(processregion(2)+(0:processregion(4)),processregion(1)+(0:processregion(3))) = crp;
  regions0 = bwlabel(regions0>0,4);
end
stat0 = regionprops(regions0,'area','boundingbox');
reg=1;
cell = 0;
regmax = length(stat0);
regold = 0;

if regmax>numel(img)/1000
  gdisp('Warning: Too many regions, consider increasing thresFactorM or threshminlevel.')
  gdisp('Some information will not be displayed ')
  regmindisp = false;
else
  regmindisp = true;
end

while reg<=regmax && reg<=p.maxRegNumber
  if reg>regold, repcount=0; else repcount=repcount+1; end
  if repcount>20, reg=reg+1; continue; end
  regold = reg;
  if regmindisp, gdisp(['processing region ' num2str(reg)]); end
  if reg>length(stat0), stat0=regionprops(regions0,'area','boundingbox'); end
  if reg>length(stat0), break; end
  statC = stat0(reg);
  
  % if the region is to small - discard it
  if statC.Area < p.areaMin
    % regions0(regions0==reg)=0; 
    if regmindisp, gdisp(['region ' num2str(reg) ' discarded, area = ' num2str(statC.Area)]); end
    reg=reg+1; 
    continue; 
  end
  
  % otherwise compute the properties for proper splitting
  roiBox(1:2) = ceil(max(statC.BoundingBox(1:2)-p.roiBorder,1)); % coordinates of the ROI box
  roiBox(3:4) = floor(min(statC.BoundingBox(1:2)+statC.BoundingBox(3:4)+p.roiBorder,[size(img,2) size(img,1)])-roiBox(1:2));
  roiRegs = imcrop(regions0,roiBox); % ROI with all regions labeled
  roiMask = bwmorph(roiRegs==reg,'close'); % ROI with region #reg labeled
  roiImg = imcrop(imge,roiBox);
  
  perim = bwperim(imdilate(roiMask,se));
  pmap = 1 - perim;
  f1=true;
  ind = 0;
  while f1 && ind<100
    ind = ind+1;
    pmap1 = 1 - perim + imerode(pmap,se);
    f1 = max(max(pmap1-pmap))>0;
    pmap = pmap1;
  end;
  regDmap = ((pmap+1) + 5*(1-roiImg/imgemax)).*roiMask;
  maxdmap = max(max(regDmap));


  % if the region is of the allowed size - try to fit the model
  if statC.Area < p.areaMax && cell<p.maxCellNumber
    
    % crop the energy and forces maps to the ROI
    roiExtDx = imcrop(extDx,roiBox);
    roiExtDy = imcrop(extDy,roiBox);
    
    % energy and forces attracting to the mask in the ROI
    % perim = bwperim(roiMask);
    % pmap = 1 - perim;
    % f1=true;
    % while f1
    %   pmap1 = 1 - perim + imerode(pmap,se);
    %   f1 = max(max(pmap1-pmap))>0;
    %   pmap = pmap1;
    % end;
    pmapEnergy = pmap + 0.1*pmap.^2;
    pmapDx = imfilter(pmapEnergy,maskdx); % distance forces
    pmapDy = imfilter(pmapEnergy,maskdy); 
    %pmapDxyMax = max(max(max(abs(pmapDx))),max(max(abs(pmapDy))));
    pmapDxyMax = 10;
    pmapEnergy = pmapEnergy/pmapDxyMax;
    pmapDx = pmapDx/pmapDxyMax; % normalize to make the max force equal to 1
    pmapDy = pmapDy/pmapDxyMax;
    
    % initial cell position
    prop = regionprops(bwlabel(roiMask),'orientation','centroid');
    theta = prop(1).Orientation*pi/180;
    x0 = prop(1).Centroid(1);
    y0 = prop(1).Centroid(2);
    
    % initial representation
    % pcCell, pcCell0 - representation of the current cell as its principal components
    if p.algorithm==2
      pcCell0 = [theta;x0;y0;zeros(p.Nkeep+1,1)];
    elseif  p.algorithm==3
      pcCell0 = [x0;y0;theta;0;zeros(p.Nkeep+1,1)];
    end

    % Making first variant of the model
    if ismember(p.algorithm,[2 3])
      % first approximation - to the exact shape of the selected region
      [pcCell,fitquality] = align(roiMask,pmapDx,pmapDy,roiExtDx*0,pcCell0,p,true,roiBox,0.5,[frame cell+1]);%0.5->thres
      gdisp(['fitquality pre-aligning = ' num2str(fitquality)])
      % adjustment of the model to the external energy map
      [pcCell,fitquality] = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,false,roiBox,thres,[frame cell+1]);
    elseif p.algorithm == 4
      pcCell = align4I(roiMask,p);
      [pcCell,fitquality] = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,roiBox,thres,[frame cell+1]);
    end
    gdisp(['fitquality aligning = ' num2str(fitquality)])
      
    % converting from box to global coordinates
    pcCellBox = pcCell;%SH130225
    pcCell = box2model(pcCell,roiBox,p.algorithm);
    % obtaining the shape of the cell in geometrical representation
    cCell = model2geom(pcCell,p.algorithm);
    
    %try splitting
    if p.algorithm==4, isct=0; else isct = intersections(cCell); end
    cellarea = statC.Area;
    if ~isempty(isct) && ~isct
      mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth);
      if length(mesh)>1
        cellarea = polyarea([mesh(:,1);flipud(mesh(:,3))],[mesh(:,2);flipud(mesh(:,4))]);
        roiMesh = mesh - repmat(roiBox([1 2 1 2])-1,size(mesh,1),1);
        res=isDivided(roiMesh,roiImg,p.splitThreshold,bgr);
        if res==-1, mesh=0; end
      end
    end
    
    % checking quality and storing the cell
    if ~isempty(pcCell) && ~isct && fitquality<p.fitqualitymax && ...
       min(cCell(:,1))>= p.imRoi(1) &&...
       min(cCell(:,2))>= p.imRoi(2) &&...
       max(cCell(:,1))<=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)) && ...
       max(cCell(:,2))<=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)) && ...
        length(mesh)>1 && p.areaMin<cellarea && p.areaMax>cellarea
      
      resplitcount = 0;

      %calculate the error estimates SH150225
      errorEstimate = getMTrackErrorEstimate(p,pcCellBox,mesh,fitquality,roiImg,roiExtDx,roiExtDy,roiMask,roiBox);
      
      % if the cell passed the quality test and it is not on the
      % boundary of the image - store it
      if p.split1 && res>0 % Splitting the cell
        mesh1 = flipud(mesh(res+1:end,:)); % daughter cell
        mesh2 = mesh(1:res-1,:); % mother cell
        
         if ismember(p.algorithm,[2 3])
          pcCell1 = splitted2model(mesh1,p);
          pcCell1 = model2box(pcCell1,roiBox,p.algorithm);
          pcCell1 = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell1,p,false,roiBox,thres,[frame cell+1]);
          pcCell1 = box2model(pcCell1,roiBox,p.algorithm);
          cCell1 = model2geom(pcCell1,p.algorithm);
          pcCell2 = splitted2model(mesh2,p);
          pcCell2 = model2box(pcCell2,roiBox,p.algorithm);
          pcCell2 = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell2,p,false,roiBox,thres,[frame cell+2]);
          pcCell2 = box2model(pcCell2,roiBox,p.algorithm);
          cCell2 = model2geom(pcCell2,p.algorithm);
        else
          pcCell1 = align4IM(mesh1,p);
          pcCell1 = model2box(pcCell1,roiBox,p.algorithm);
          cCell1 = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell1,p,roiBox,thres,[frame cell+1]);
          cCell1 = box2model(cCell1,roiBox,p.algorithm); % Corrected pcCell->cCell 2008/08/02
          pcCell2 = align4IM(mesh2,p);
          pcCell2 = model2box(pcCell2,roiBox,p.algorithm);
          cCell2 = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell2,p,roiBox,thres,[frame cell+2]);
          cCell2 = box2model(cCell2,roiBox,p.algorithm); % Corrected pcCell->cCell 2008/08/02
        end

        mesh1 = model2mesh(cCell1,p.meshStep,p.meshTolerance,p.meshWidth);
        mesh2 = model2mesh(cCell2,p.meshStep,p.meshTolerance,p.meshWidth);
        cellStruct1.algorithm = p.algorithm;
        cellStruct2.algorithm = p.algorithm;
        cellStruct1.birthframe = frame;
        cellStruct2.birthframe = frame;
        if size(pcCell1,2)==1, model=cCell1'; else model=cCell1; end
        cellStruct1.model = model;
        if size(pcCell2,2)==1, model=cCell2'; else model=cCell2; end
        cellStruct2.model = model;
        cellStruct1.mesh = mesh1;
        cellStruct2.mesh = mesh2;
        cellStruct1.polarity = 1;
        cellStruct2.polarity = 1;
        cellStruct1.stage = 1;
        cellStruct2.stage = 1;
        cellStruct1.timelapse = tl;
        cellStruct2.timelapse = tl;
        cellStruct1.divisions = [];
        cellStruct2.divisions = [];
        cellStruct1.box = roiBox;
        cellStruct2.box = roiBox;
        cellStruct1.ancestors = [];
        cellStruct2.ancestors = [];
        cellStruct1.descendants = [];
        cellStruct2.descendants = [];
        cellStruct1.errorEstimate = errorEstimate;
        cellStruct2.errorEstimate = errorEstimate;
        cell = cell+2;
        cellList{frame}{cell-1} = cellStruct1;
        cellList{frame}{cell} = cellStruct2;
        gdisp(['fitting region ' num2str(reg) ' - passed and saved as cells ' num2str(cell-1) ' and ' num2str(cell)])
        reg=reg+1;
        continue;
      end
      
      % if the cell passed the quality test and it is not on the
      % boundary of the image - store it
      cellStruct.algorithm = p.algorithm;
      cellStruct.birthframe = frame;
      if size(pcCell,2)==1, model=pcCell'; else model=pcCell; end
      cellStruct.model = model;
      cellStruct.mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth);
      cellStruct.polarity = 0;
      cellStruct.ancestors = [];
      cellStruct.descendants = [];
      cellStruct.stage = 1;
        %cellStruct.stage = getStage(roiImg,cellStruc.model);
      cellStruct.timelapse = 0;
      cellStruct.divisions = [];
      cellStruct.box = roiBox;
      cellStruct.errorEstimate = errorEstimate;
      cell = cell+1;
      cellList{frame}{cell} = cellStruct;
      gdisp(['fitting region ' num2str(reg) ' - passed and saved as cell ' num2str(cell)])
      reg=reg+1;
      continue;
    end
    % if the cell is not on the image border OR it did not pass the
    % quality test - split it
    reason = 'unknown';
    if isempty(pcCell), reason = 'no cell found'; 
    elseif fitquality>=p.fitqualitymax, reason = 'bad fit quality'; 
    elseif min(cCell(:,1))<=p.imRoi(1), reason = 'cell on x=0 boundary';
    elseif min(cCell(:,2))<=p.imRoi(2), reason = 'cell on y=0 boundary';
    elseif max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)), reason = 'cell on x=max boundary';
    elseif max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)), reason = 'cell on y=max boundary';
    elseif isct, reason = 'model has intersections';
    elseif length(mesh)<=1, reason = 'problem getting mesh'; 
    elseif p.areaMin>=cellarea, reason = 'cell too small'; 
    elseif p.areaMax<=cellarea, reason = 'cell too big'; 
    end
    if isempty(who('resplitcount')), resplitcount=0; end
    if (p.areaMin>cellarea) || (resplitcount>=10) % Discarding cells
      gdisp(['fitting region ' num2str(reg) ' - quality check failed - ' reason])
      reg=reg+1;
      resplitcount = 0;
      continue;
    elseif isfield(p,'splitbndcells') && ~p.splitbndcells && ...
         ( min(cCell(:,1))<=p.imRoi(1) || min(cCell(:,2))<=p.imRoi(2) || ... 
            max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)) || ...
            max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)))
      gdisp(['fitting region ' num2str(reg) ' - ' reason ' - splitting not allowed'])
      reg=reg+1;
      resplitcount = 0;
      continue;
    else
      resplitcount = resplitcount+1;
      gdisp(['fitting region ' num2str(reg) ' - quality check failed - try resplitting (' num2str(resplitcount) ')'])
    end
  end
  
  if ~p.splitregions % do not split, just discard
    reg=reg+1;
    gdisp(['region ' num2str(reg) ' discarded - splitting not allowed'])
    continue;
  end
  
  % If the area is too large or if the adjustment failed - splitting the image

%  regWshed = roiMask.*watershed(max(max(regDmap))-regDmap,4);
%   roiLabeled = roiMask;
%   for i=1:maxdmap
%     %roiLabeled = bwlabel(imerode(roiLabeled,se));
%     roiLabeled = roiLabeled.*(regDmap>i);
%     nsreg = max(max(roiLabeled));
%     if nsreg<=1, continue; end
%     if nsreg>=2, break; end
%   end
%   [yL,xL]=ind2sub(roiBox([4 3])+1,find(roiLabeled==1,1));
%   roiLabeled = ~((regWshed==0).*(roiLabeled~=1));
%   roiLabeled = bwselect(roiLabeled,xL,yL,4)==1;
%   roiLabeled = imdilate(roiLabeled,se).*roiMask;
%   
%   statL=regionprops(bwlabel(roiLabeled,8),'area');
%   roiResidual = roiMask.*~roiLabeled;
%   statM=regionprops(bwlabel(roiResidual,8),'area');

  [roiLabeled,roiResidual] = splitonereg(roiMask.*roiImg);
  if isempty(roiLabeled), reg=reg+1; gdisp(['region ' num2str(reg) ' discarded - unable to split (1)']); continue; end
  statL=regionprops(bwlabel(roiLabeled),'area');
  statM=regionprops(bwlabel(roiResidual),'area');
  
  if isempty(statL), area1 = -1; 
  elseif statL(1).Area<p.areaMin, area1 = 0;
  elseif statL(1).Area<p.areaMax, area1 = 1;
  else area1 = 2;
  end

  if isempty(statM), area2 = -1;
  elseif statM(1).Area<p.areaMin, area2 = 0;
  elseif statM(1).Area<p.areaMax, area2 = 1;
  else area2 = 2;
  end
  
%   if i==maxdmap, area1=0; area2=0; end

  if area1==-1 || area2==-1
    % roiRegs = roiRegs.*(~roiMask);
    reg=reg+1;
    gdisp(['region ' num2str(reg) ' discarded - unable to split (2)'])
    continue;
  end  
  
  if area1>=1 && area2>=1
    regmax = regmax+1;
    roiRegs = roiRegs.*(~roiMask) + reg*roiLabeled + regmax*roiResidual;
    roiMask = roiLabeled;
    gdisp(['region ' num2str(reg) ' splitted and the number of regions increased - return for test'])
  elseif area1==0 && area2>=1
    roiRegs = roiRegs.*(~roiMask) + reg*roiResidual;
    roiMask = roiResidual;
    area1 = area2;
  elseif area2==0
    roiRegs = roiRegs.*(~roiMask) + reg*roiLabeled;
    roiMask = roiLabeled;
    gdisp(['region ' num2str(reg) ' splitted - return for test'])
  end

  if area1==0
    roiRegs = roiRegs.*(~roiMask);
    reg=reg+1;
    gdisp(['region ' num2str(reg) ' discarded'])
  end
  if area1>0
    statF = regionprops(bwlabel(roiRegs==reg),'Area');
    stat0(reg).Area=statF(1).Area;
  end
  
  regions0(roiBox(2):roiBox(2)+roiBox(4),roiBox(1):roiBox(1)+roiBox(3)) = roiRegs;
  continue;
end
cellListN(frame) = length(cellList{frame});
end % processFrameI function













function proccells2 = processFrameOld(frame,proccells,cont,processregion)
% This function is for processing the first frame only
% 
% It searches for the cells on the frame and then does refinement

global imsizes se cellList cellListN rawPhaseData p;

% In ~cont mode the output array proccells2 stays empty.
% In the cont mode, proccells2 will stay empty if the input array is empty
% (but all cells will be processed), otherwise it will hold the progeny of 
% the proccells. In the cont mode, regardless of the proccells, the progeny
% is determined by positioning of the center of the new cell inside of the
% outline of the old cell. The first cell matched this way is mother, the
% second is daughter, the rest are discarded.
proccells2 = [];
if ~isempty(proccells), lst = proccells; elseif frame>1, lst = 1:length(cellList{frame-1}); else lst = []; end

if p.invertimage
  if isempty(rawPhaseData)||frame>size(rawPhaseData,3), gdisp(['The frame ' num2str(frame) ' will not be processed: no image']); return; end
  img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
else
  if isempty(rawPhaseData)||frame>size(rawPhaseData,3), gdisp(['The frame ' num2str(frame) ' will not be processed: no image']); return; end
  img = rawPhaseData(:,:,frame);
end

imge = img2imge(img,p.erodeNum);
imge16 = img2imge16(img,p.erodeNum);
imgemax = max(max(imge));
img = img2imge(img,0);

thres = graythreshreg(imge,p.threshminlevel);
regions0 = getRegions(p.edgemode,imge,thres,p.thresFactorM,p.opennum,imge16,p.edgeSigmaL,p.edgeSigmaV,p.valleythresh1,p.valleythresh2,p.imRoi);
if ~isempty(processregion)
  crp = regions0(processregion(2)+(0:processregion(4)),processregion(1)+(0:processregion(3)));
  regions0 = regions0*0;
  regions0(processregion(2)+(0:processregion(4)),processregion(1)+(0:processregion(3))) = crp;
  regions0 = bwlabel(regions0>0,4);
end
stat0 = regionprops(regions0,'area','boundingbox');
reg=1;
cell = 0;
regmax = length(stat0);
regold = 0;

if regmax>numel(img)/500
  gdisp('Too many regions, consider increasing thresFactorM.')
  gdisp('Some information will not be displayed ')
  regmindisp = false;
else
  regmindisp = true;
end

while reg<=regmax && reg<=p.maxRegNumber
  if reg>regold, repcount=0; else repcount=repcount+1; end
  if repcount>20, reg=reg+1; continue; end
  regold = reg;
  if regmindisp, gdisp(['processing region ' num2str(reg)]); end
  if reg>length(stat0), stat0=regionprops(regions0,'area','boundingbox'); end
  if reg>length(stat0), break; end
  statC = stat0(reg);
  
  % if the region is to small - discard it
  if statC.Area < p.areaMin
    % regions0(regions0==reg)=0; 
    if regmindisp, gdisp(['region ' num2str(reg) ' discarded, area = ' num2str(statC.Area)]); end
    reg=reg+1; 
    continue; 
  end
  
  % otherwise compute the properties for proper splitting
  roiBox(1:2) = max(statC.BoundingBox(1:2)-p.roiBorder,1); % coordinates of the ROI box
  roiBox(3:4) = min(statC.BoundingBox(1:2)+statC.BoundingBox(3:4)+p.roiBorder,[size(img,2) size(img,1)])-roiBox(1:2);
  roiRegs = imcrop(regions0,roiBox); % ROI with all regions labeled
  roiMask = roiRegs==reg; % ROI with region #reg labeled
  % for k=1:p.dilMaskNum, roiMask = imdilate(roiMask,se); end
  roiImg = imcrop(imge,roiBox);
  roiImgI = imcrop(img,roiBox);
  
  perim = bwperim(imdilate(roiMask,se));
  pmap = 1 - perim;
  f1=true;
  while f1
    pmap1 = 1 - perim + imerode(pmap,se);
    f1 = max(max(pmap1-pmap))>0;
    pmap = pmap1;
  end;
  regDmap = (3*(pmap+1) + 5*(1-roiImg/imgemax)).*roiMask;
  maxdmap = max(max(regDmap));


  % if the region is of the allowed size - accept it
  if statC.Area < p.areaMax
    
    % Get and smooth the boundary of the region
    [ii,jj]=find(bwperim(roiMask),1,'first');
    cCell=bwtraceboundary(roiMask,[ii,jj],'n',8,inf,'counterclockwise');
    if isfield(p,'interpoutline') && p.interpoutline
      ctrlist = interpoutline(cCell,roiImg,p);
    else
      ctrlist = {cCell};
    end
    if isempty(ctrlist), gdisp(['fitting region ' num2str(reg) ' - failed, no outline creating during smoothing']); end
    for cind = 1:length(ctrlist)
      cCell = ctrlist{cind};
      if p.fsmooth>0 && p.fsmooth<Inf
        fpp = frdescp(cCell);
        cCell = ifdescp(fpp,p.fsmooth);
        cCell=[cCell;cCell(1,:)]; % make the last point the same as the first one
      end
      cCell = rot90(cCell,2);
      cCell(:,1) = cCell(:,1)+roiBox(1)-1;
      cCell(:,2) = cCell(:,2)+roiBox(2)-1;
      
      % Save the data
      if p.getmesh, mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth); end
      if (~p.getmesh || length(mesh)>4) && ...
         min(cCell(:,1))<=p.imRoi(1) && min(cCell(:,2))<=p.imRoi(2) && ...
         max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)) && ...
         max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1))
        cellStruct.algorithm = 1;
        if p.getmesh, cellStruct.mesh = mesh; end
        cellStruct.polarity = 0;
        cellStruct.box = roiBox;
        cellStruct.ancestors = [];
        cellStruct.descendants = [];
        cellStruct.divisions = [];
        cellStruct.stage = [];
        cellStruct.contour = cCell;
        cellStruct.polarity = 0; % algorithm 1 cannot divide cells unless this is interpoutline
        if cont && frame>1
          if p.getmesh
            xcent = sum(mesh(floor(size(mesh,1)/2),[1 3]))/2;
            ycent = sum(mesh(floor(size(mesh,1)/2),[2 4]))/2;
          else
            xcent = mean(cCell(:,1));
            ycent = mean(cCell(:,2));
          end
          mother = [];
          for mcell=lst % length(cellList{frame})<mcell ||
            if isempty(cellList{frame-1}{mcell}), continue; end
            if ~isfield(cellList{frame-1}{mcell},'contour') && (~isfield(cellList{frame-1}{mcell},'mesh') ...
                || length(cellList{frame-1}{mcell}.mesh)<=4), continue; end
            if length(cellList{frame-1}{mcell}.box)>1
              box = cellList{frame-1}{mcell}.box;
              if inpolygon(xcent,ycent,[box(1) box(1) box(1)+box(3) box(1)+box(3)],[box(2) box(2)+box(4) box(2)+box(4) box(2)])
                if isfield(cellList{frame-1}{mcell},'mesh')
                  mesh2 = cellList{frame-1}{mcell}.mesh;
                  ctrx = [mesh2(:,1);flipud(mesh2(:,3))];
                  ctry = [mesh2(:,2);flipud(mesh2(:,4))];
                elseif isfield(cellList{frame-1}{mcell},'contour')
                  ctrx = cellList{frame-1}{mcell}.contour(:,1);
                  ctry = cellList{frame-1}{mcell}.contour(:,2);
                end
                if inpolygon(xcent,ycent,ctrx,ctry)
                  mother=mcell;
                  break;
                end
              end
            end
          end
          cell = cell+1;
          if isempty(mother)
            gdisp(['fitting region ' num2str(reg) ' - no parent found - discarded'])
          else
            cellStruct.birthframe = cellList{frame-1}{mcell}.birthframe;
            if length(cellList{frame})<mother || isempty(cellList{frame}{mother})
              cellList{frame}{mother} = cellStruct;
              if cont && ~isempty(proccells),  proccells2 = [proccells2 mother]; end
              gdisp(['fitting region ' num2str(reg) ' - passed and saved'])
            else
              cellList{frame}{mother}.divisions = [cellList{frame}{mother}.divisions frame];
              daughter = getdaughter(mother,length(cellList{frame}{mother}.divisions),cellListN(frame-1));
              cellList{frame}{mother}.descendants = daughter;
              if length(cellList{frame})<daughter || isempty(cellList{frame}{daughter})
                if daughter<p.maxCellNumber
                  cellList{frame}{daughter} = cellStruct;
                  if cont && ~isempty(proccells), proccells2 = [proccells2 daughter]; end
                else
                  gdisp(['fitting region ' num2str(reg) ' - new cell still born - overpopulation'])
                end
              else
                gdisp(['fitting region ' num2str(reg) ' - too many children, one discarded'])
              end
            end
          end
        else
          cellStruct.birthframe = frame;
          cell = cell+1;
          cellList{frame}{cell} = cellStruct;
          gdisp(['fitting region ' num2str(reg) ' - passed and saved as cell ' num2str(cell)])
        end
      else
        % if the cell is not on the image border OR it did not pass the
        % quality test - split it
        reason = 'unknown';
        if p.getmesh && length(mesh)<=4, reason = 'bad region shape quality'; end
        if min(cCell(:,1))<=p.imRoi(1), reason = 'cell on x=0 boundary'; end
        if min(cCell(:,2))<=p.imRoi(2), reason = 'cell on y=0 boundary'; end
        if max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)), reason = 'cell on x=max boundary'; end
        if max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)), reason = 'cell on y=max boundary'; end
        gdisp(['fitting region ' num2str(reg) ' - quality check failed - ' reason])
      end
    end
    reg=reg+1;
    continue
  elseif ~p.splitregions % do not split, just discard
    reg=reg+1;
    gdisp(['region ' num2str(reg) ' discarded - splitting not allowed'])
    continue
  end
  
  % If the area is too large or if the adjustment failed - splitting the image

%  regWshed = roiMask.*watershed(max(max(regDmap))-regDmap,4);
%   roiLabeled = roiMask;
%   for i=1:maxdmap
%     %roiLabeled = bwlabel(imerode(roiLabeled,se));
%     roiLabeled = roiLabeled.*(regDmap>i);
%     nsreg = max(max(roiLabeled));
%     if nsreg<=1, continue; end
%     if nsreg>=2, break; end
%   end
%   [yL,xL]=ind2sub(roiBox([4 3])+1,find(roiLabeled==1,1));
%   roiLabeled = ~((regWshed==0).*(roiLabeled~=1));
%   roiLabeled = bwselect(roiLabeled,xL,yL,4)==1;
%   roiLabeled = imdilate(roiLabeled,se).*roiMask;
%   
%   statL=regionprops(bwlabel(roiLabeled,8),'area');
%   roiResidual = roiMask.*~roiLabeled;
%   statM=regionprops(bwlabel(roiResidual,8),'area');

  [roiLabeled,roiResidual] = splitonereg(roiMask.*roiImgI);
  if isempty(roiLabeled), reg=reg+1; gdisp(['region ' num2str(reg) ' discarded - unable to split (1)']); continue; end
  statL=regionprops(roiLabeled ,'area');
  statM=regionprops(roiResidual,'area');

  if isempty(statL), area1 = -1; 
  elseif statL(1).Area<p.areaMin, area1 = 0;
  elseif statL(1).Area<p.areaMax, area1 = 1;
  else area1 = 2;
  end

  if isempty(statM), area2 = -1;
  elseif statM(1).Area<p.areaMin, area2 = 0;
  elseif statM(1).Area<p.areaMax, area2 = 1;
  else area2 = 2;
  end
  
  % if i==maxdmap, area1=0; area2=0; end

  if area1==-1 || area2==-1
    reg=reg+1;
    gdisp(['region ' num2str(reg) ' discarded - unable to split'])
    continue;
  end  
  
  if area1>=1 && area2>=1
    regmax = regmax+1;
    roiRegs = roiRegs.*(~roiMask) + reg*roiLabeled + regmax*roiResidual;
    roiMask = roiLabeled;
    gdisp(['region ' num2str(reg) ' splitted and the number of regions increased - return for test'])
  elseif area1==0 && area2>=1
    roiRegs = roiRegs.*(~roiMask) + reg*roiResidual;
    roiMask = roiResidual;
    area1 = area2;
  elseif area2==0
    roiRegs = roiRegs.*(~roiMask) + reg*roiLabeled;
    roiMask = roiLabeled;
    gdisp(['region ' num2str(reg) ' splitted - return for test'])
  end

  if area1==0
    roiRegs = roiRegs.*(~roiMask);
    reg=reg+1;
    gdisp(['region ' num2str(reg) ' discarded'])
  end
  if area1>0
    statF = regionprops(bwlabel(roiRegs==reg),'Area');
    stat0(reg).Area=statF(1).Area;
  end
  
  regions0(roiBox(2):roiBox(2)+roiBox(4),roiBox(1):roiBox(1)+roiBox(3)) = roiRegs;
end
if cont
  cellListN(frame) = cellListN(frame-1);
else
  cellListN(frame) = length(cellList{frame});
end
end

% --- End of major processing functions ---















% --- Processing subfunctions ---

function [pcc,ftq] = align(im,mpx,mpy,ngmap,pcc,p,f,roiBox,thres,celldata)
% A modification of align function with adding non-linear bending
% works with p.algorithm == 2 or 3
% 
% Parameters: 
% p.fitDisplay1, p.fitConvLevel1, p.fitMaxIter1, p.fitStep1 - 1st frame
% p.fitDisplay, p.fitConvLevel, p.fitMaxIter, p.fitStep - refine step
% f = true - region fitting (frame=1), else - refinement step (frame>=1)
% 
% New model array (pcc) includes:
% 1   - rotation
% 2   - bending
% 3,4 - shift
% 5   - stretching along x
% 6-  - other parameters deternimed by the principal components analysis

  global coefPCA coefS N weights dMax;
  pcc = reshape(pcc,[],1);
  if f % fitting to "mask"
    dsp = p.fitDisplay1;
    lev = p.fitConvLevel1;
    amax = p.fitMaxIter1;
    stp = p.fitStep1;
  else % fitting to the force field
    dsp = p.fitDisplay;
    lev = p.fitConvLevel;
    amax = p.fitMaxIter;
    stp = p.fitStep;
  end
  if dsp
    fig = createdispfigure([celldata f]);
    nextstop = 1;
    contmode = false;
  else
    nextstop = Inf;
  end
  
  weights2 = weights.^p.rigidity;
  Kstp = 1;
      % if p.algorithm==2, Kstp=ones(4+p.Nkeep,1); end
      % if p.algorithm==3, Kstp=ones(5+p.Nkeep,1); end

  [cCell,cCell2] = model2geom(pcc,p.algorithm);
  if isempty(cCell), pcc=[]; ftq=0; return; end
  for a=1:amax
    % Compute forces
    xCell = cCell(:,1);
    yCell = cCell(:,2);
    Fx = interp2a(1:roiBox(3)+1,1:roiBox(4)+1,mpx,xCell,yCell,'linear',0);
    Fy = interp2a(1:roiBox(3)+1,1:roiBox(4)+1,mpy,xCell,yCell,'linear',0);

    if ~f % for alignment to the image
      Tx = -(circshift(cCell(:,2),-1) - circshift(cCell(:,2),1));
      Ty = -(circshift(cCell(:,1),1) - circshift(cCell(:,1),-1));
      dtxy = sqrt(Tx.^2+Ty.^2);%sqrt(sum(Tx.^2 + Ty.^2)/N);
      Tx = Tx./dtxy;
      Ty = Ty./dtxy;
      TxM = repmat(cCell(:,1),1,p.attrRegion+1)+Tx*(0:p.attrRegion);
      TyM = repmat(cCell(:,2),1,p.attrRegion+1)+Ty*(0:p.attrRegion);
      Clr0 = interp2a(1:roiBox(3)+1,1:roiBox(4)+1,im,TxM,TyM,'linear',0);
      Clr = 1-1./(1+(Clr0/thres/p.thresFactorF).^p.attrPower);
      are = polyarea(cCell(:,1),cCell(:,2));
      Tnr = - p.neighRepA * interp2a(1:roiBox(3)+1,1:roiBox(4)+1,ngmap,xCell,yCell,'linear',0);
      % T = p.attrCoeff * sum(Clr,2) - p.repCoeff * (are<p.repArea*p.areaMax) * (1-Clr(:,1)) + Tnr;
      T = p.attrCoeff * mean(Clr(:,p.attrRegion+1:end),2) - p.repCoeff * (are<p.repArea*p.areaMax) * (1-mean(Clr(:,1:p.attrRegion+1),2)) + Tnr;
      Tx = Tx.*T;
      Ty = Ty.*T;
      F = [-Fx;Fy] + [Tx;Ty]; % F is a vector of forces in the image frame of reference
    else % for pre-alignment to the mask
      Tx = -(circshift(cCell(:,2),-1) - circshift(cCell(:,2),1));
      Ty = -(circshift(cCell(:,1),1) - circshift(cCell(:,1),-1));
      dtxy = sqrt(sum(Tx.^2 + Ty.^2)/N);
      Tx = Tx/dtxy;
      Ty = Ty/dtxy;
      TxM = repmat(cCell(:,1),1,p.attrRegion+1)+Tx*(0:p.attrRegion);
      TyM = repmat(cCell(:,2),1,p.attrRegion+1)+Ty*(0:p.attrRegion);
      Clr0 = interp2a(1:roiBox(3)+1,1:roiBox(4)+1,im,TxM,TyM,'linear',0);
      Clr = 1-1./(1+(Clr0/thres/p.thresFactorF).^p.attrPower); % ?
      are = polyarea(cCell(:,1),cCell(:,2));
      T = p.attrCoeff1 * sum(Clr,2) - p.repCoeff1 * (are<p.repArea*p.areaMax) * (1-Clr(:,1));
      Tx = Tx.*T;
      Ty = Ty.*T;
      F = [-Fx;Fy] + [Tx;Ty]; % F is a vector of forces in the image frame of reference
    end
    if p.algorithm==2
      F2 = M(reshape(F,[],2),-pcc(1));
      F2 = [F2(:,1);F2(:,2)];
      Fpc = weights.*([coefPCA(:,1:2)'*F;coefPCA(:,3:end)'*F2]);

      % rigidity forces
      Fpc = Fpc - p.rigidity*[0;0;0;ones(p.Nkeep,1)].*(pcc(2:end).^3).*(1./(weights.^2)-1);

      %Fpc = ([coefPCA(:,1:2)'*F;coefPCA(:,3:end)'*F2]);
      if max(Fpc.^2)>0, Fpc = Fpc/max(Fpc.^2).^0.25; end

      % Compute torque
      Mm = sum( - Fx.*(yCell-mean(yCell)) - Fy.*(xCell-mean(xCell)) )/N/dMax/dMax;
      pcc(1) = pcc(1) + Mm;
      
      if a>1, Fpc2old=Fpc2; end
      Fpc2 = [Mm;Fpc];
    elseif p.algorithm==3
      % Check 1
      F2 = M(reshape(F,[],2),-pcc(3)); % F2 is a vector of forces in the cell frame

      F3 = Bf(cCell2,F2,pcc(4)); % F3 is a vector of forces in the unbent cell frame
      F3y = F3(:,2); % Only y componet affects the bending in unbent cell frame

      F3 = [F3(:,1);F3(:,2)];
      Fpc = weights(5:end).*(coefPCA'*F3);
      % Fpc is a vector of forces acting on the linear components
      
      % if f % Fitting to mask ???????????
      %   Fpc(1) = Fpc(1) + 0*weights(5);
      % end
      
      % Add rigidity forces
      % Fpc = Fpc - p.rigidity*[0;ones(p.Nkeep,1)].*(pcc(5:end).^3).*(1./(weights(5:end).^2)-1);
      pcc = pcc.*weights2; % Imposing rigidity, modified from a force method 08/18/08

      % Compute torque
      Mm = sum( - Fx.*(yCell-mean(yCell)) - Fy.*(xCell-mean(xCell)) )/N/dMax/dMax;

      % Compute bending "torque", must be done in the unbent cell frame
      Bm = -sum((F3y-mean(F3y)).*cCell2(:,1).^2)/N/dMax/dMax/dMax;

      % Normalize Fpc
      if a>1, Fpc2old=Fpc2; end
      Fpc2 = [weights(1:4).*[coefS'*F;Mm;Bm];Fpc];
      if max(Fpc.^2)>0, Fpc2 = Fpc2/max(Fpc2.^2).^0.25; end
    end
    
    if a>1
      K = sum(sign(Fpc2.*Fpc2old))/(5+p.Nkeep);
      if K<0.3
        Kstp=Kstp/1.5;
      elseif K>0.7
        Kstp=min(1,Kstp.*1.5);
      end
    end
    % Shift & transform
    pcc = pcc + Kstp*stp*p.scaleFactor*Fpc2;
    [cCell,cCell2] = model2geom(pcc,p.algorithm);
    
    % if a>1
    %   K = sign(Fpc2.*Fpc2old);
    %   Kstp(K>0) = min(1,Kstp(K>0).*1.5);
    %   Kstp(K<0) = Kstp(K<0)/1.5;
    % end
    % Shift & transform
    % pcc = pcc + Kstp.*stp.*Fpc2;
    % [cCell,cCell2] = model2geom(pcc,p.algorithm);
    
    if a>=nextstop % Displaying results
      figure(fig)
      %emap = emap-min(min(emap));
       % image(repmat(im,[1 1 3])/max(max(im+0)));
      %colormap(gray(2));
      imshow(im,[]);
      set(gca,'nextplot','add');
      plot(cCell(:,1),cCell(:,2));
      for b=1:N/2-1; plot(cCell([b N-b],1),cCell([b N-b],2),'m');end
      quiver(cCell(:,1),cCell(:,2),-3*Fx,3*Fy,0,'r');
      %quiver(cCell(:,1),cCell(:,2),Tx,Ty,0,'y');
      set(gca,'nextplot','replace');
      setdispfiguretitle(fig,celldata,a)
      drawnow
      if ~contmode || a==p.fitMaxIter
        waitfor(fig,'UserData');
        if ishandle(fig)
          u = get(fig,'UserData');
        else
          u = 'stop';
        end
        if strcmp(u,'next')
          nextstop = a+1;
        elseif strcmp(u,'next100')
          nextstop = a+100;
        elseif strcmp(u,'skip')
          nextstop = Inf;
        elseif strcmp(u,'continue')
          nextstop = a+1;
          contmode = 1;
        elseif strcmp(u,'stop')
          % ftq = sum(abs(Fpc2))/length(Fpc2);
          % break;
          error('Testing mode terminated')
        elseif strcmp(u,'debug')
           dbstop if warning MicrobeTracker:DebugStop
           warning('MicrobeTracker:DebugStop','Program stopped in debugger')
        end
        set(fig,'UserData','');
      end
      pause(0.05);
    end

    % Condition to finish
    ftq = sum(abs(Fpc2))/length(Fpc2);
    if ftq<lev; break; end
  end
  pcc = reshape(pcc,1,[]);
  % evaluate the result by some "closeness" criterium
  % fitquality - is a weighter mean square deviation from the average
  % cell in pixels. Typical values: 1-2
  % TODO: The rejection value should be estimated from image analysis
  % ftq =sum((interp2a(1:roiBox(3)+1,1:roiBox(4)+1,emap,xCell,yCell,'linear',0)).^2);
  % ftq = sum((interp2a(1:roiBox(3)+1,1:roiBox(4)+1,emap,xCell,yCell,'linear',0)).^4);
end

function cCell = align4I(msk,p)
  % Get and smooth the boundary of the region
  [ii,jj]=find(bwperim(msk),1,'first');
  pp=bwtraceboundary(msk,[ii,jj],'n',4,inf,'counterclockwise');
  fpp = frdescp(pp);
  cCell = ifdescp(fpp,p.fsmooth);
  mesh = model2mesh(cCell,p.fmeshstep,p.meshTolerance,p.meshWidth);
  if length(mesh)>4
    cCell = fliplr([mesh(:,1:2);flipud(mesh(2:end-1,3:4))]);
  else
    cCell = [];
  end
  cCell = makeccw(cCell);
end

function b = makeccw(a)
  if isempty(a)
    b = [];
  else
    if isContourClockwise(a)
      b = circshift(flipud(a),1);
    else
      b=a;
    end
  end
end

function fig = createdispfigure(celldata)
  global dhandles
  if exist('dhandles','var') && isfield(dhandles,'fig') && ishandle(dhandles.fig)
    fig = dhandles.fig;
    setdispfiguretitle(dhandles.fig,celldata,0)
    set(dhandles.fig,'UserData','');
    return
  end 
  dhandles.fig = figure('KeyPressFcn',@dispfitkeypress,'CloseRequestFcn',@dispfitclosereq,'Toolbar','none','Menubar','none','NumberTitle','off','IntegerHandle','off','UserData','');
  setdispfiguretitle(dhandles.fig,celldata,0)
  pos = get(dhandles.fig,'pos');
  pos = [pos(1)+pos(3)/2-200 pos(2)+pos(4)/2-300 400 430];
  set(dhandles.fig,'pos',pos);
  fig = dhandles.fig;
  dhandles.ax = axes('units','pixels','pos',[1 1 400 400],'box','off','tick','out','DataAspectRatio',[1 1 1]);
  dhandles.cpanel = uipanel(dhandles.fig,'units','pixels','pos',[1 401 400 30]);
  dhandles.next = uicontrol(dhandles.cpanel,'units','pixels','Position',[5 4 63 20],'String','Next step','callback',@dispfitcontrol,'KeyPressFcn',@dispfitkeypress);
  dhandles.next100 = uicontrol(dhandles.cpanel,'units','pixels','Position',[70 4 63 20],'String','+100 steps','callback',@dispfitcontrol,'KeyPressFcn',@dispfitkeypress);
  dhandles.skip = uicontrol(dhandles.cpanel,'units','pixels','Position',[135 4 63 20],'String','Skip cell','callback',@dispfitcontrol,'KeyPressFcn',@dispfitkeypress);
  dhandles.continue = uicontrol(dhandles.cpanel,'units','pixels','Position',[200 4 63 20],'String','Continue','callback',@dispfitcontrol,'KeyPressFcn',@dispfitkeypress);
  dhandles.debug = uicontrol(dhandles.cpanel,'units','pixels','Position',[265 4 63 20],'String','Debug','callback',@dispfitcontrol,'KeyPressFcn',@dispfitkeypress);
  dhandles.stop = uicontrol(dhandles.cpanel,'units','pixels','Position',[330 4 63 20],'String','Stop','callback',@dispfitcontrol,'KeyPressFcn',@dispfitkeypress);
  set(dhandles.fig,'ResizeFcn',@dispfitresizefcn);
  dhandles.ax;
  
  function dispfitkeypress(hObject,eventdata)
    if isempty(eventdata) || ~isfield(eventdata,'Key')
      return
    elseif strcmp(eventdata.Key,'rightarrow')
      set(dhandles.fig,'UserData','next');
    elseif strcmp(eventdata.Key,'downarrow')
      set(dhandles.fig,'UserData','next100');
    elseif strcmp(eventdata.Key,'escape')
      set(dhandles.fig,'UserData','stop');
    end
  end
  function dispfitcontrol(hObject,eventdata)
    if hObject==dhandles.next
      set(dhandles.fig,'UserData','next');
    elseif hObject==dhandles.next100
      set(dhandles.fig,'UserData','next100');
    elseif hObject==dhandles.skip
      set(dhandles.fig,'UserData','skip');
    elseif hObject==dhandles.continue
      set(dhandles.fig,'UserData','continue');
    elseif hObject==dhandles.stop
      set(dhandles.fig,'UserData','stop');
    elseif hObject==dhandles.debug
      set(dhandles.fig,'UserData','debug');
    end
  end
  function dispfitclosereq(hObject, eventdata)
    delete(dhandles.fig)
    clear('dhandles')
  end
  function dispfitresizefcn(hObject, eventdata)
    screenSize = get(0,'ScreenSize');
    pos = get(dhandles.fig,'pos');
    pos = [pos(1:2) max(pos(3:4),[400 430])];
    pos(2) = min(pos(2),screenSize(4)-460);
    set(dhandles.fig,'pos',pos);
    set(dhandles.ax,'pos',[1 1 pos(3) pos(4)-30]);
    set(dhandles.cpanel,'pos',[1 pos(4)-30+1 pos(3) 30]);
  end
end

function setdispfiguretitle(fig,celldata,a)
  if ~isempty(celldata)
    if celldata(1)==0, frame='?'; else frame = num2str(celldata(1)); end
    if celldata(2)==0, cell='?'; else cell = num2str(celldata(2)); end
    if length(celldata)<3 || celldata(3)==0, mode='fit'; else mode = 'mask'; end
    nm = ['Aligning: frame ' frame ', cell ' cell ', step ' num2str(a) ', ' mode ' mode'];
  else
    nm = ['Aligning: step: ' num2str(a)];
  end
  set(fig,'Name',nm);
end

function cCell = align4IM(mesh,p)
  % Align to a 'bad' mesh
%   pp=[mesh(:,1:2);flipud(mesh(2:end-1,3:4))];
%   fpp = frdescp(pp);
%   cCell = ifdescp(fpp,p.fsmooth);
%   mesh = model2mesh(cCell,p.fmeshstep);
  if length(mesh)>4
    cCell = [mesh(:,1:2);flipud(mesh(2:end-1,3:4))];
  else
    cCell = [];
  end
  cCell = makeccw(cCell);
end
  
function [pcc2,ftq] = align4(im,mpx,mpy,ngmap,pcc,p,roiBox,thres,celldata)
% A modification of align function to work with filamentous cells
% works only with p.algorithm == 4, f is absolete
% Parameters: 
% p.fitDisplay, p.fitConvLevel, p.fitMaxIter, p.fitStep - refine step
% Algorithm 4 model array (pcc) is 2L array

  % Simplifying some p parameters
  if ~isfield(p,'maxmesh'), p.maxmesh = 1000; end
  isnanpcc = isnan(sum(pcc,2));
  if max(isnanpcc)==1, pcc=pcc(~isnanpcc,:); end
  if isempty(pcc) || length(pcc)>p.maxmesh*4
    pcc2 = [];
    ftq = 0;
    return
  end
  dsp = p.fitDisplay;
  if dsp
    fig = createdispfigure(celldata);
    nextstop = 1;
    contmode = false;
  else
    nextstop = Inf;
  end

  % Get and smooth the boundary of the region
  if ~(isfield(p,'smoothbeforefit') && ~p.smoothbeforefit)
    fpp = frdescp(pcc);
    cCell = ifdescp(fpp,p.fsmooth);
    mesh = model2mesh(cCell,p.fmeshstep,p.meshTolerance,p.meshWidth);
    if length(mesh)>4
      pcc = [mesh(:,1:2);flipud(mesh(2:end-1,3:4))];
    else
      pcc = [];
    end
    pcc = makeccw(pcc);
  end
  if isempty(pcc) || isempty(im)
    pcc2 = [];
    ftq = 0;
    return
  end
  
  % Initializations
  L = size(pcc,1);
  N = L/2+1;
  stp = 1;
  H = ceil(p.cellwidth*pi/2/stp/2);
  %heads = [ones(1,H),zeros(1,N-2*H),ones(1,2*H-1),zeros(1,N-2*H),ones(1,H-1)]';
  %hcorr = heads * ddx*pi/4/(H-1);
  xCell = pcc(:,1);
  yCell = pcc(:,2);
  Kstp = 1;
  
  % Construct model cell and the rigidity forces in it
  ddx = round(p.rigidityRange);
  A = 1/2./(1:ddx);
  A = A/(2*sum(A)-A(1));
  
  ddxB = round(p.rigidityRangeB);
  B = 1/2./sqrt(1:ddxB);
  B = B/(2*sum(B)-B(1));
  
  HA = H+1:-1:1;
  HA = pi*HA/(2*sum(HA)-HA(1)); % ???
  
  x(H+1+2*ddx) = 0;
  y(H+1+2*ddx) = 0;
  for i=H+2*ddx:-1:H+1
    x(i) = x(i+1) - stp;
    y(i) = y(i+1);
  end
  alpha = HA(H+1);
  for i=H:-1:1
    x(i) = x(i+1) - stp*cos(alpha);
    y(i) = y(i+1) + stp*sin(alpha);
    alpha = HA(i) + alpha;
  end
  x = [fliplr(x(2:end)),x];
  y = [2*y(1)-fliplr(y(2:end)),y];
  y = y*p.cellwidth*p.scaleFactor/abs(y(1));
  [fx,fy] = getrigidityforces(x',y',A);
  [Tx,Ty] = getnormals(x',y');
  f = Tx.*fx + Ty.*fy;
  f = f((end+1)/2:(end+1)/2+H);
  hcorr = [f;zeros(N-2*H-2,1);flipud(f);f(2:end);zeros(N-2*H-2,1);flipud(f(2:end))];
  rgt = [1:H+1 (H+1)*ones(1,N-2*H-2) H+1:-1:1 2:H+1 (H+1)*ones(1,N-2*H-2) H+1:-1:2]'/(H+1);
  if length(rgt)>L % N-2
    L2 = length(rgt); 
    rm = (L2-L)/2;
    lv = ceil(N/2)-1;
    hcorr = hcorr([1:lv lv+1+rm:L2/2+1+lv L2/2+1+lv+1+rm:end]);
    rgt = rgt([1:lv lv+1+rm:L2/2+1+lv L2/2+1+lv+1+rm:end]); 
  end
  
  % Opposite points interaction (Fdstx,Fdsty)
  lside = 2:N-1;
  rside = L:-1:N+1;
  % wdt = [zeros(1,H) H*ones(1,N-2*H-2) zeros(1,H)]'/H;
  cellWidth = [2*abs(y(end/2+3/2:end/2+1/2+H)-y(end/2+1/2))';p.cellwidth*p.scaleFactor*ones(N-2*H-2,1);2*abs(y(end/2+1/2+H:-1:end/2+3/2)-y(end/2+1/2))'];
  if length(cellWidth)>N-2, cellWidth = cellWidth([1:ceil(N/2)-1 end+2-floor(N/2):end]); end
  wdt = cellWidth/max(cellWidth);
  
  for a=1:p.fitMaxIter
    % Vector image forces (Fix,Fiy)
    Fix =-p.imageforce * interp2a(1:roiBox(3)+1,1:roiBox(4)+1,mpx,xCell,yCell,'linear',0);
    Fiy = p.imageforce * interp2a(1:roiBox(3)+1,1:roiBox(4)+1,mpy,xCell,yCell,'linear',0);
    
    % Get normals to the centerline (Tx,Ty)
    Tx = circshift(yCell,-1) - circshift(yCell,1);
    Ty = circshift(xCell,1) - circshift(xCell,-1);
    dtxy = sqrt(Tx.^2 + Ty.^2);
    dtxy = dtxy + mean(dtxy)/100;
    if min(dtxy)==0, pcc2=[];ftq=0;return;end
    Tx = Tx./dtxy;
    Ty = Ty./dtxy;
    % [Tx,Ty] = getnormals(xCell,yCell);
    
    Lx = Ty;
    Ly = -Tx;
    
    % Get area outside the cell & attraction/repulsion (Fax,Fay)
      % Matrix attRegion wide outside
    TxM = repmat(xCell,1,2*p.attrRegion+1)+Tx*(-p.attrRegion:p.attrRegion);
    TyM = repmat(yCell,1,2*p.attrRegion+1)+Ty*(-p.attrRegion:p.attrRegion);
      % Matrix attRegion wide outside
    Clr0 = interp2a(1:roiBox(3)+1,1:roiBox(4)+1,im,TxM,TyM,'linear',0);
      % Non-normalized 'area' attraction
    Clr = 1-1./(1+(Clr0/thres/p.thresFactorF).^p.attrPower);
      % Cell area
    are = polyarea(xCell,yCell);
      % Scalar repulsion/attraction forces
    T = p.attrCoeff * mean(Clr(:,p.attrRegion+1:end),2) - p.repCoeff * (are<p.repArea*p.areaMax) * (1-mean(Clr(:,1:p.attrRegion+1),2));
      % Vector repulsion/attraction forces
    Fax = Tx.*T;
    Fay = Ty.*T;
    
    T = - p.neighRepA * interp2a(1:roiBox(3)+1,1:roiBox(4)+1,ngmap,xCell,yCell,'linear',0);
    Fnrx = Tx.*T;
    Fnry = Ty.*T;
    
    % Opposite points interaction (Fdstx,Fdsty)
    Dst = sqrt((xCell(lside)-xCell(rside)).^2 + (yCell(lside)-yCell(rside)).^2);
    Fdst = (cellWidth-Dst)./cellWidth;
    g = 5;
    Fdst((Dst./cellWidth)<0.5)=Fdst((Dst./cellWidth)<0.5).*g-(g-1)*0.5;
    Fdst = p.wspringconst*wdt.*Fdst.*cellWidth;
    Fdst1x = Fdst.*(xCell(lside)-xCell(rside))./Dst;
    Fdst1y = Fdst.*(yCell(lside)-yCell(rside))./Dst;
    Fdstx = zeros(L,1);
    Fdsty = zeros(L,1);
    Fdstx(lside) = Fdst1x;
    Fdsty(lside) = Fdst1y;
    Fdstx(rside) = -Fdst1x;
    Fdsty(rside) = -Fdst1y;
    
    % Rigidity (Frx,Fry)
    [D4x,D4y] = getrigidityforces(xCell,yCell,A);
    Frx = p.rigidity * (D4x - Tx.*hcorr).*rgt;
    Fry = p.rigidity * (D4y - Ty.*hcorr).*rgt;
    
    % Backbone rigidity (Fbrx,Fbry)
    xCnt = (xCell(1:N)+flipud(xCell([N:L 1])))/2;
    yCnt = (yCell(1:N)+flipud(yCell([N:L 1])))/2;
    Tbx = circshift(yCnt,-1) - circshift(yCnt,1);
    Tby = circshift(xCnt,1) - circshift(xCnt,-1);
    Tbx(end) = Tbx(end-1);
    Tby(end) = Tby(end-1);
    Tbx(1) = Tbx(2);
    Tby(1) = Tby(2);
    dtxy = sqrt(Tbx.^2 + Tby.^2);
    Tbx = Tbx./dtxy;
    Tby = Tby./dtxy;
    [D4btx,D4bty] = getrigidityforcesL(xCnt,yCnt,B);
    D4b = (D4btx.*Tbx + D4bty.*Tby)/2;
    D4bx = p.rigidityB * D4b.*Tbx;
    D4by = p.rigidityB * D4b.*Tby;
    Fbrx = [D4bx;flipud(D4bx(2:end-1))];
    Fbry = [D4by;flipud(D4by(2:end-1))];

    % Perpendicular ribs (Fpx,Fpy)
    Fpbx = xCell(lside)-xCell(rside);
    Fpby = yCell(lside)-yCell(rside);
    Fp = Fpbx.*Tbx(2:end-1) + Fpby.*Tby(2:end-1);
    Fpx = zeros(L,1);
    Fpy = zeros(L,1);
    Fpbx = p.horalign*(Fpbx-Fp.*Tbx(2:end-1));
    Fpby = p.horalign*(Fpby-Fp.*Tby(2:end-1));
    Fpx(lside) = -Fpbx;
    Fpy(lside) = -Fpby;
    Fpx(rside) = Fpbx;
    Fpy(rside) = Fpby;
    
    % Equal distances between points (Fqx,Fqy)
    Fqx = p.eqaldist*(circshift(xCell,1)+circshift(xCell,-1)-2*xCell);
    Fqy = p.eqaldist*(circshift(yCell,1)+circshift(yCell,-1)-2*yCell);
    Fq = Lx.*Fqx + Ly.*Fqy;
    Fqx = Fq.*Lx;
    Fqy = Fq.*Ly;
    
    % Get the resulting force
    if a>1, Fo = [Fx;Fy]; end
    if isempty(who('Kstp2'))||isempty(Kstp2), Kstp2 = 1; end
    Fx = (Fix + Fax + Fnrx + Fdstx) + Frx;
    Fy = (Fiy + Fay + Fnry + Fdsty) + Fry;
    Fs = Fx.*Tx + Fy.*Ty;
    Fx = Fs.*Tx + Fpx + Fbrx + Fqx;
    Fy = Fs.*Ty + Fpy + Fbry + Fqy;

    % Normalize
    Fm = abs(Fs).^0.2;
    Fm = Fm+mean(Fm)/100;
    if min(Fm)==0, gdisp('Error - zero force'); pcc2=[]; ftq = 0; return; end

    if a>1
      K = sum((Fo.*[Fx;Fy])>0)/2/L;
      if K<0.4
        Kstp=Kstp/1.4; % Kstp saturates oscillations perpendicular to the coutour
      elseif K>0.6
        Kstp=min(1,Kstp.*1.2);
      end
    end
    mxf = p.fitStep*Kstp*max(max(abs(Fx))+max(abs(Fy)));
    asd = (xCell-circshift(xCell,1)).^2+(yCell-circshift(yCell,1)).^2;
    mnd = sqrt(min(asd));
    med = sqrt(mean(asd));
    Kstp2 = min([Kstp2*1.1 1 mnd/mxf/2 3*mnd/med]); % Katp2 prevents points crossing
    if isfield(p,'moveall') && p.moveall>0
      if a>1
        mfxold = mfx;
        mfyold = mfy;
        MFold = MF;
      end
      xm = xCell-mean(xCell);
      ym = yCell-mean(yCell);
      MF = mean(-Fy.*xm+Fx.*ym);
      MI = sum(xm.^2+ym.^2);
      Fmrx =  ym*MF/MI;
      Fmry = -xm*MF/MI;
      mfx = mean(Fx);
      mfy = mean(Fy);
      if isfield(p,'fitStepM'), mfitstep = p.fitStepM; else mfitstep = p.fitStep; end
      if isempty(who('Kstpm'))||isempty(Kstpm), Kstpm = mfitstep; end
      Kstpm = min([Kstpm*1.5 mfitstep/(sqrt(mean(Fx)^2+mean(Fy)^2)) mfitstep/abs(MF)*sqrt(MI)]); % Katpm prevents large (>1px) mean steps
      if a>1 && (mfx*sign(mfxold)<-abs(mfxold)/2 || mfy*sign(mfyold)<-abs(mfyold)/2 || MF*sign(MFold)<-abs(MFold)/2)
        Kstpm = Kstpm/2;
      end
    end
    
    if a>=nextstop % Displaying results
      figure(fig)
      imshow(im,[]);
      set(gca,'nextplot','add');
      plot(xCell,yCell,'k');
      plot(xCnt,yCnt,'k');
      for b=2:L/2; plot(xCell([b L+2-b]),yCell([b L+2-b]),'m');end
      %quiver(xCell,yCell,3*Fx,3*Fy,0,'r');
      %quiver(xCell,yCell,3*Tx,3*Ty,0,'color',[0 0 0]);
      %quiver(xCnt,yCnt,3*Tbx,3*Tby,0,'color',[0 0 0]);
      quiver(xCell,yCell,3*Fix,3*Fiy,0,'color',[1 0 0]);
      quiver(xCell,yCell,3*Fax,3*Fay,0,'color',[0.7 0 0.7]);
      quiver(xCell,yCell,3*Frx,3*Fry,0,'color',[0 0 1]);
      quiver(xCell,yCell,3*Fdstx,3*Fdsty,0,'color',[1 0.5 0]);
      %quiver(xCell,yCell,3*p.horalign*(Fdstx-Fp.*Tx),3*p.horalign*(Fdsty-Fp.*Ty),0,'color',[1 0.3 .3]);
      quiver(xCell,yCell,3*Fqx,3*Fqy,0,'color',[0.6 0.9 0]);
      quiver(xCnt,yCnt,3*D4bx,3*D4by,0,'color',[0 0.9 0.6]);
      quiver(xCell,yCell,3*Fpx,3*Fpy,0,'color',[1 0.9 0.2]);
      quiver(xCell,yCell,3*Fx,3*Fy,0,'color',[0.3 0.3 1]);
      %quiver(cCell(:,1),cCell(:,2),Tx,Ty,0,'y');
      set(gca,'nextplot','replace');
      setdispfiguretitle(fig,celldata,a)
      drawnow
      if ~contmode || a==p.fitMaxIter
        waitfor(fig,'UserData');
        if ishandle(fig)
          u = get(fig,'UserData');
        else
          u = 'stop';
        end
        if strcmp(u,'next')
          nextstop = a+1;
        elseif strcmp(u,'next100')
          nextstop = a+100;
        elseif strcmp(u,'skip')
          nextstop = Inf;
        elseif strcmp(u,'continue')
          nextstop = a+1;
          contmode = 1;
        elseif strcmp(u,'stop')
          % ftq = sum(abs(Fm))/length(Fm)/2;
          % break;
          error('Testing mode terminated')
        elseif strcmp(u,'debug')
           dbstop if warning MicrobeTracker:DebugStop
           warning('MicrobeTracker:DebugStop','Program stopped in debugger')
        end
        set(fig,'UserData','');
      end
      pause(0.05);
    end
    
    % Move
    if isfield(p,'moveall') && p.moveall>0
      Fx = Kstp*Kstp2*p.scaleFactor*p.fitStep*Fx*(1-p.moveall)+Kstpm*(mean(Fx)+Fmrx)*p.moveall;
      Fy = Kstp*Kstp2*p.scaleFactor*p.fitStep*Fy*(1-p.moveall)+Kstpm*(mean(Fy)+Fmry)*p.moveall;
    else
      Fx = Kstp*Kstp2*p.scaleFactor*p.fitStep*Fx;
      Fy = Kstp*Kstp2*p.scaleFactor*p.fitStep*Fy;
    end
    xCell = xCell + Fx;
    yCell = yCell + Fy;
    
    if max(isnan(xCell))==1
      pcc2=[];
      ftq=0;
      return
    end
    
    % Looking for self-intersections
    [i1,i2]=intxySelfC(xCell,yCell);
    % Moving points halfway to the projection on the opposite strand
    iMovCurveArr = [];
    xMovCurveArr = [];
    yMovCurveArr = [];
    for i=1:2:(length(i1)-1)
      if i1(i)<=i1(i+1)
        iMovCurve = mod((i1(i)+1:i1(i+1))-1,L)+1;
      else
        iMovCurve = mod((i1(i)+1:i1(i+1)+L)-1,L)+1;
      end
      if length(iMovCurve)<2, continue; end
      if i2(i)+1>=i2(i+1)
        iRefCurve = mod((i2(i)+1:-1:i2(i+1))-1,L)+1;
      else
        iRefCurve = mod((i2(i)+1+L:-1:i2(i+1))-1,L)+1;
      end
      % iMovCurve = mod((i1(i)+1:i1(i+1))-1,L)+1;
      % if length(iMovCurve)<2, continue; end
      % iRefCurve = mod((i2(i)+1:-1:i2(i+1))-1,L)+1;
      xMovCurve = reshape(xCell(iMovCurve),1,[]);
      yMovCurve = reshape(yCell(iMovCurve),1,[]);
      xRefCurve = reshape(xCell(iRefCurve),1,[]);
      yRefCurve = reshape(yCell(iRefCurve),1,[]);
      [xMovCurve,yMovCurve]=projectCurve(xMovCurve,yMovCurve,xRefCurve,yRefCurve);
      iMovCurveArr = [iMovCurveArr iMovCurve];
      xMovCurveArr = [xMovCurveArr xMovCurve];
      yMovCurveArr = [yMovCurveArr yMovCurve];
    end
    xCell(iMovCurveArr) = xMovCurveArr;
    yCell(iMovCurveArr) = yMovCurveArr;
    
    % Condition to finish
    ftq = sum(abs(Fm))/length(Fm)/2;
    if ftq < p.fitConvLevel; break; end
  end

  % Output model
  pcc2 = [xCell,yCell];
  
  if isfield(p,'smoothafterfit') && p.smoothafterfit
    fpp = frdescp(pcc2);
    cCell = ifdescp(fpp,p.fsmooth);
    mesh = model2mesh(cCell,p.fmeshstep,p.meshTolerance,p.meshWidth);
    if length(mesh)>4
      pcc2 = [mesh(:,1:2);flipud(mesh(2:end-1,3:4))];
    else
      pcc2 = [];
    end
    pcc2 = makeccw(pcc2);
  end
end

function [fx,fy] = getrigidityforces(x,y,A)
  fx = - 2*sum(A)*x;
  fy = - 2*sum(A)*y;
  for i=1:length(A)
    fx = fx + A(i)*(circshift(x,i)+circshift(x,-i));
    fy = fy + A(i)*(circshift(y,i)+circshift(y,-i));
  end
end

function [fx,fy] = getrigidityforcesL(x,y,A)
  %A=0.5;
  L=length(A);
  fx = 0*x;
  fy = 0*y;
%   l = length(x);
%   D = sqrt((x([1 2 l-2 l-1])-x([2 3 l-1 l])).^2+(y([1 2 l-2 l-1])-y([2 3 l-1 l])).^2);
%   x(1) = x(2) + (x(1)-x(2))*D(2)/D(1);
%   y(1) = y(2) + (y(1)-y(2))*D(2)/D(1);
%   x(end) = x(end-1) + (x(end)-x(end-1))*D(3)/D(4);
%   y(end) = y(end-1) + (y(end)-y(end-1))*D(3)/D(4);
  for i=1:L
    fxt = A(i)*(x(1:end-2*i)/2+x(2*i+1:end)/2-x(i+1:end-i));
    fyt = A(i)*(y(1:end-2*i)/2+y(2*i+1:end)/2-y(i+1:end-i));
    fx(i+1:end-i) = fx(i+1:end-i) + fxt;
    fy(i+1:end-i) = fy(i+1:end-i) + fyt;
    fx(1:end-2*i) = fx(1:end-2*i) - fxt/2;
    fy(1:end-2*i) = fy(1:end-2*i) - fyt/2;
    fx(2*i+1:end) = fx(2*i+1:end) - fxt/2;
    fy(2*i+1:end) = fy(2*i+1:end) - fyt/2;
  end
end

function [tx,ty] = getnormals(x,y)
  tx = circshift(y,-1) - circshift(y,1);
  ty = circshift(x,1) - circshift(x,-1);
  dtxy = sqrt(tx.^2 + ty.^2);
  tx = tx./dtxy;
  ty = ty./dtxy;
end








% *********** Conversion functions, conputational subfunctions ************

function pcc = model2box(pcc,box,alg)
  % converts model ("pcc") from global to local (determind by "box")
  % coordinates, algorithm ("alg")-specific
   if length(pcc)<=1, return; end
   if alg==2
     pcc(2:3) = pcc(2:3) - box(1:2) + 1;
   elseif alg==3
     pcc(1:2) = pcc(1:2) - box(1:2) + 1;
   elseif alg==4
     pcc(:,1) = pcc(:,1) - box(1) + 1;
     pcc(:,2) = pcc(:,2) - box(2) + 1;
   end
end

function pcc = box2model(pcc,box,alg)
  % converts model ("pcc") from local to global (determind by "box")
  % coordinates, algorithm ("alg")-specific
   if isempty(pcc), return; end
   if alg==2
     pcc(2:3) = pcc(2:3) + box(1:2) - 1;
   elseif alg==3
     pcc(1:2) = pcc(1:2) + box(1:2) - 1;
   elseif alg==4
     pcc(:,1) = pcc(:,1) + box(1) - 1;
     pcc(:,2) = pcc(:,2) + box(2) - 1;
   end
end

function [cell,cell2] = model2geom(pcc,alg)
% cell - cell geometry
% cell2 - unbent & unrotated cell geometry
global coefPCA mCell
% This function switches paramertic representation of the model to
% geometrical one
  if length(pcc)<=1, cell=[];cell2=[]; return; end
  if (alg==2 || alg==3) && (isempty(coefPCA) || isempty(mCell)), cell=[];cell2=[]; return; end
  if alg==2
    % pcc structure: [theta, x, y, size, Nkeep parameters from PCA]
    cell2 = mCell+reshape(coefPCA*[0;0;reshape(pcc(4:end),[],1)],[],2);
    cell = M(cell2,pcc(1));
    cell(:,1) = cell(:,1)+pcc(2);
    cell(:,2) = cell(:,2)+pcc(3);
  elseif alg==3
    % pcc structure: [x, y, theta, curvature, length, Nkeep parameters from PCA]
    cell2 = mCell+reshape(coefPCA*reshape(pcc(5:end),[],1),[],2); % ???
    cellT = B(cell2,pcc(4)); % cellT = B(cell2,pcc(4));
    cell = M(cellT,pcc(3));
    cell(:,1) = cell(:,1)+pcc(1);
    cell(:,2) = cell(:,2)+pcc(2);
  elseif alg==4
    cell = pcc;
    cell2 = pcc;
  end
end

function b=B(a,t)
  % bends a set of points (down if t>0) with the radius of curvature 1/t
  if abs(t)<1/10000000, b=a; return; end
  R = sign(t)*max(1/abs(t),1.1*max(-a(:,2)));
  b(:,1)=(a(:,2)+R).*sin(a(:,1)/R);
  b(:,2)=(a(:,2)+R).*cos(a(:,1)/R)-R;
end

function f2=Bf(c,f,t)
  % converts forces from a bend to non-bent frame
  % now c and f2 a are in the non-bent configuration and f is in the bent one
  % essentially the point rotates by angle phi
  if abs(t)<1/10000000, f2=f; return; end
  R = sign(t)*max(1/abs(t),1.1*max(-c(:,2)));
  phi = c(:,1)/R;
  cs = cos(phi);
  sn = sin(phi);
  
  f2(:,1) = cs.*f(:,1) - sn.*f(:,2);
  f2(:,2) = sn.*f(:,1) + cs.*f(:,2);
end

function a=M(a,t)
  % rotates a set of points clockwise by an angle t and scales it by a factor s
  cost = cos(t);
  sint = sin(t);
  mt = [cost -sint; sint cost];
  for i=1:size(a,1)
    a(i,:)=a(i,:)*mt;
  end
end



function res=intersections(coord)
  % determines if a contour (Nx2 array of x-y pairs of coordinates) is
  % self-intersecting
  if isempty(coord), res = []; return; end
  alpha = 0.01;
  c1 = (coord+circshift(coord,1))/2;
  c2 = alpha*(coord-circshift(coord,1));
  IN = inpolygon(-c1(:,1).*c2(:,2),c1(:,2).*c2(:,1),coord(:,1),coord(:,2));
  if length(IN)~=sum(IN) && sum(IN)~=0, res=true; else res=false; end
  %if find(IN,1), res=true; else res=false; end
end



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
    imRoi = varargin{6};
end

if ~exist('imRoi','var')
   imRoi = [1 1 Inf Inf];
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
imgProc = cropToRoi(imgProc, imRoi); %SH 130226
res = bwlabel(imgProc,4);
end


function [extDx,extDy]=getExtForces(im,im16,p)
% computes the image forces from the image im using parameters
% p.forceWeights, p.dmapThres, p.dmapPower, and p.gradSmoothArea
% "im" - original phase image, "im16" - eroded image, outpits energy mast -
% same size as the microscope images.

global maskdx maskdy
try % TODO: make sure try-catch is necessary, also profile the speed
  
  % Edge forces
  [fex,fey] = edgeforce(im16,p.edgeSigmaL);

  % Gradient forces
  gradSmoothFilter = fspecial('gaussian',2*ceil(1.5*p.gradSmoothArea)+1,p.gradSmoothArea);
  gradSmoothImage = imfilter(im,gradSmoothFilter); % this filtered image is going to be used in ther sections
  gradThres = mean(mean(gradSmoothImage(gradSmoothImage>mean(mean(gradSmoothImage)))));
  gradEnergy = imfilter(im,maskdx).^2 + imfilter(im,maskdy).^2;
  gradEnergy = - gradEnergy./(gradEnergy + gradThres^2);
  gradDx = imfilter(gradEnergy,maskdx);
  gradDy = imfilter(gradEnergy,maskdy);
  gradDxyMax = max(max(max(abs(gradDx))),max(max(abs(gradDy))));
  gradEnergy = gradEnergy/gradDxyMax; % normalize to make the max force equal to 1

  % Threshold forces
  thresLevel = p.thresFactorF*graythreshreg(gradSmoothImage,p.threshminlevel); % Use the same filtering as in gradient section
  thresEnergy = (gradSmoothImage-thresLevel).^2;
  thresDx = imfilter(thresEnergy,maskdx);
  thresDy = imfilter(thresEnergy,maskdy);
  thresDxyMax = max(max(max(abs(thresDx))),max(max(abs(thresDy))));
  thresEnergy = thresEnergy/thresDxyMax; % normalize to make the max force equal to 1

  % Combined force
  extEnergy = gradEnergy*p.forceWeights(2) ...
        + thresEnergy*p.forceWeights(3);
  extDx = imfilter(extEnergy,maskdx)+fex*p.forceWeights(1);
  extDy = imfilter(extEnergy,maskdy)+fey*p.forceWeights(1);
catch
  extDx=[];
  extDy=[];
  le = lasterror;
  gdisp('Error in getExtForces function')
  gdisp(['Error message: ' le.message '. Error in line ' le.stack.line])
  gdisp('Most probable reasons of the error - memory leaks or missing parameters. Try restarting MATLAB.')
end
end



function im2 = img2imge16(im,nm)
% erodes image "im" by "nm" pixels
global se
im2 = im;
for i=1:nm, im2 = imerode(im2,se);end
end



function im2 = img2imge(im,nm)
% erodes image "im" by "nm" pixels and normalizes it, outputs double
global se
im2 = im;
for i=1:nm, im2 = imerode(im2,se);end
im2 = double(im2);
mn=mmin(im2);
mx=mmax(im2);
im2=1-(im2-mn)/double(mx-mn);
end



function res=isDivided(mesh,img,thr,bgr)
% splits the cell based on phase profile
% takes the mesh, cropped phase image (img) and threshold - the minimum
% relative drop in the profile to consider the cell divided
% 
% Current problem: the profile is based on the colors already in the phase
% image, which makes the output dependent of the contrast of the image
if length(mesh)<5, res=0; return; end
isz = size(img);
lng = sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2);%1-img/max(max(img))
% box = [1 1 size(img,2) size(img,1)];
% prf = getOneSignal(mesh,box,img,1);

img = max(img-bgr,0);

prf = interp2a(1:isz(2),1:isz(1),img,(mesh(:,1)+mesh(:,3))/2,(mesh(:,2)+mesh(:,4))/2,'linear',0);
prf = 0.5*prf + 0.25*(prf([1 1:end-1])+prf([2:end end]));
prf = prf.*lng;

mn = mean(prf);
minima = [false reshape( (prf(2:end-1)<prf(1:end-2))&(prf(2:end-1)<=prf(3:end))|(prf(2:end-1)<=prf(1:end-2))&(prf(2:end-1)<prf(3:end)) ,1,[]) false];
if isempty(minima) || sum(prf)==0, res=-1; return; end
while true
  if sum(minima)==0, res=0; return; end
  im = find(minima);
  [min0,i0] = min(prf(minima));
  max1 = max(prf(1:im(i0)));
  max2 = max(prf(im(i0):end));
  if max1<0.5*mn || max2<0.5*mn, minima(im(i0))=0; continue; else break; end
end
if thr>=1
  res = 0;
elseif (max1+max2-2*min0)/(max1+max2)>thr
  res = im(i0);
else
  res = 0;
end
end

function joincells(frame,lst)
  % this function looks for the cells close enough according to
  % parameters and tries to join all such pairs
  
  global cellList p se rawPhaseData maskdx maskdy

  if checkparam(p,'invertimage','algorithm','erodeNum','joindist','joinangle','roiBorder','meshStep','meshTolerance','meshWidth','splitThreshold')
    gdisp('Joining cells failed: one or more required parameters not provided.');
    return
  end
  if frame>size(rawPhaseData,3), gdisp('Joining cells failed: unacceptable frame number'); return; end
  if p.invertimage
    img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
  else
    img = rawPhaseData(:,:,frame);
  end
  imge = img2imge(img,p.erodeNum);
  sz = size(img);
  % p.joindist = 5;
  % p.joinangle = 0.2;
  if length(cellList)<frame || isempty(cellList{frame}), return; end
  if isempty(lst), lst=1:length(cellList{frame}); end
  joindist = p.joindist^2;
  
  imge16 = img2imge16(img,p.erodeNum);
  thres = graythreshreg(imge,p.threshminlevel);
  bgr = phasebgr(imge,thres);
  [extDx,extDy] = getExtForces(imge,imge16,p);
  if isempty(extDx), gdisp('Joining cells failed: unable to get energy'); return; end
  
  for cell1=lst
    if isempty(cellList{frame}{cell1}), continue; end
    mesh1 = cellList{frame}{cell1}.mesh;
    if length(mesh1)<=4, continue; end
    for ori1 = [1 size(mesh1,1)]
      if isempty(cellList{frame}{cell1}), continue; end
      for cell2=lst
        if cell2<cell1
          if isempty(cellList{frame}{cell2}), continue; end
          mesh2 = cellList{frame}{cell2}.mesh;
          if length(mesh2)<=4, continue; end
          for ori2=[1 size(mesh2,1)]
            d2 = (mesh1(ori1,1)-mesh2(ori2,1))^2+(mesh1(ori1,2)-mesh2(ori2,2))^2;
            if d2<joindist
              a1 = angle(2*mesh1(ori1,1)-mesh1(abs(ori1-3)+1,1)-mesh1(abs(ori1-3)+1,3) + 1i*(2*mesh1(ori1,2)-mesh1(abs(ori1-3)+1,2)-mesh1(abs(ori1-3)+1,4)));
              a2 = angle(2*mesh2(ori2,1)-mesh2(abs(ori2-3)+1,1)-mesh2(abs(ori2-3)+1,3) + 1i*(2*mesh2(ori2,2)-mesh2(abs(ori2-3)+1,2)-mesh2(abs(ori2-3)+1,4)));
              a = pi/2-abs(mod(a1-a2,pi)-pi/2);
              if abs(a)<p.joinangle
                % Now actuallty do the joining procedire
                border = p.roiBorder;
                % roiBox = round([min(min([mesh1(:,[1 3]);mesh2(:,[1 3])]))-border min(min([mesh1(:,[2 4]);mesh2(:,[2 4])]))-border...
                %         max(max([mesh1(:,[1 3]);mesh2(:,[1 3])]))-min(min([mesh1(:,[1 3]);mesh2(:,[1 3])]))+2*border...
                %         max(max([mesh1(:,[2 4]);mesh2(:,[2 4])]))-min(min([mesh1(:,[2 4]);mesh2(:,[2 4])]))+2*border]);
                roiBox = [max(round(min(min([mesh1(:,[1 3]);mesh2(:,[1 3])]))-border),1) max(round(min(min([mesh1(:,[2 4]);mesh2(:,[2 4])]))-border),1)];
                roiBox(3) = min(round(max(max([mesh1(:,[1 3]);mesh2(:,[1 3])]))+border),sz(2)-1)-roiBox(1);
                roiBox(4) = min(round(max(max([mesh1(:,[2 4]);mesh2(:,[2 4])]))+border),sz(1)-1)-roiBox(2);
                mask1 = poly2mask([mesh1(:,1);flipud(mesh1(:,3))]-roiBox(1)+1,[mesh1(:,2);flipud(mesh1(:,4))]-roiBox(2)+1,roiBox(4)+1,roiBox(3)+1);
                mask2 = poly2mask([mesh2(:,1);flipud(mesh2(:,3))]-roiBox(1)+1,[mesh2(:,2);flipud(mesh2(:,4))]-roiBox(2)+1,roiBox(4)+1,roiBox(3)+1);
                mask3 = poly2mask([mesh1(abs(ori1-2)+1,1) mesh2(abs(ori2-2)+1,1) mesh2(abs(ori2-2)+1,3) mesh1(abs(ori1-2)+1,3)]-roiBox(1)+1,...
                  [mesh1(abs(ori1-3)+1,2) mesh2(abs(ori2-3)+1,2) mesh2(abs(ori2-3)+1,4) mesh1(abs(ori1-3)+1,4)]-roiBox(2)+1,roiBox(4)+1,roiBox(3)+1);
                mask4 = poly2mask([mesh1(abs(ori1-2)+1,1) mesh2(abs(ori2-2)+1,3) mesh2(abs(ori2-2)+1,1) mesh1(abs(ori1-2)+1,3)]-roiBox(1)+1,...
                  [mesh1(abs(ori1-3)+1,2) mesh2(abs(ori2-3)+1,4) mesh2(abs(ori2-3)+1,2) mesh1(abs(ori1-3)+1,4)]-roiBox(2)+1,roiBox(4)+1,roiBox(3)+1);
                mask = imdilate(min(1,mask1+mask2+mask3+mask4),se);
                edg = bwperim(mask);

                pmap = 1 - edg;
                f1=true;
                while f1
                  pmap1 = 1 - edg + imerode(pmap,se);
                  f1 = max(max(pmap1-pmap))>0;
                  pmap = pmap1;
                end;
                pmapEnergy = pmap + 0.1*pmap.^2;
                pmapDx = imfilter(pmapEnergy,maskdx); % distance forces
                pmapDy = imfilter(pmapEnergy,maskdy); 
                pmapDxyMax = 10;
                pmapDx = pmapDx/pmapDxyMax; % normalize to make the max force equal to 1
                pmapDy = pmapDy/pmapDxyMax;
                
                % roiBox2 = roiBox+[1 1 0 0];%[roiBox(2) roiBox(1) roiBox(4)-1 roiBox(3)-1]; % standard format
                roiImg = imcrop(imge,roiBox);
                roiExtDx = imcrop(extDx,roiBox);
                roiExtDy = imcrop(extDy,roiBox);
                roiBox(3:4) = [size(roiExtDx,2) size(roiExtDx,1)]-1; % !
                
                prop = regionprops(bwlabel(mask),'orientation','centroid');
                theta = prop(1).Orientation*pi/180;
                x0 = prop(1).Centroid(1);
                y0 = prop(1).Centroid(2);

                if p.algorithm==2
                  pcCell0 = [theta;x0;y0;zeros(p.Nkeep+1,1)];
                elseif p.algorithm==3
                  pcCell0 = [x0;y0;theta;0;zeros(p.Nkeep+1,1)];
                end

                if ismember(p.algorithm,[2 3])
                  [pcCell,fitquality] = align(mask,pmapDx,pmapDy,pmapDx*0,pcCell0,p,true,roiBox,0.5,[frame cell2]);
                  [pcCell,fitquality] = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,false,roiBox,thres,[frame cell2]);
                elseif p.algorithm == 4
                  pcCell = align4I(mask,p);
                  [pcCell,fitquality] = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,roiBox,thres,[frame cell2]);
                end
                pcCell = box2model(pcCell,roiBox,p.algorithm);
                % if p.algorithm==2
                %   pcCell(2) = pcCell(2) + roiBox(1);
                %   pcCell(3) = pcCell(3) + roiBox(2);
                % elseif  p.algorithm==3
                %   pcCell(1) = pcCell(1) + roiBox(1);
                %   pcCell(2) = pcCell(2) + roiBox(2);
                % end
                cCell = model2geom(pcCell,p.algorithm);
                mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth);
                
                res=isDivided(mesh-repmat(roiBox(1:2),size(mesh,1),2),roiImg,p.splitThreshold,bgr);
                if length(mesh)<5 || res || isempty(res)
                  gdisp(['Tried to join cells ' num2str(cell1) ' and ' num2str(cell2) ' - failed, resplit again'])
                  continue
                else
                  cellList{frame}{cell2}.mesh = mesh;
                  if size(pcCell,2)==1, model=reshape(pcCell,[],1); else model=pcCell; end
                  cellList{frame}{cell2}.model = model; % TODO: check if works for alg. 2-3
                  cellList{frame}{cell2}.box = roiBox;
                  cellList{frame}{cell1} = [];
                  gdisp(['Tried to join cells ' num2str(cell1) ' and ' num2str(cell2) ' - success, cell ' num2str(cell1) ' removed'])
                end
              end
            end
          end
        end
      end
    end
  end    
end

function refinecell(frame,lst)
  % this function runs the alignment again for the selected cells
  % only for algorithms 2-4
  
  global cellList cellListN p rawPhaseData
  
  if checkparam(p,'invertimage','algorithm','erodeNum','meshStep','meshTolerance','meshWidth')
    gdisp('Refining cells failed: one or more required parameters not provided.');
    return
  end
  if length(cellList)<frame || isempty(cellList{frame}) || isempty(lst), return; end
  if ~ismember(p.algorithm,[2 3 4]), gdisp('There is no refinement routine for algorithm 1'); return; end
  if frame>size(rawPhaseData,3), return; end
  if p.invertimage    
    img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
  else
    img = rawPhaseData(:,:,frame);
  end
  imge = img2imge(img,p.erodeNum);
  imge16 = img2imge16(img,p.erodeNum);
  thres = graythreshreg(imge,p.threshminlevel);
  [extDx,extDy] = getExtForces(imge,imge16,p);
  if isempty(extDx), gdisp('Refining cells failed: unable to get energy'); return; end
  n = 0;
  for cell=lst
    prevStruct = cellList{frame}{cell};
    roiBox = prevStruct.box;
    roiImg = imcrop(imge,roiBox);
    roiBox(3:4) = [size(roiImg,2) size(roiImg,1)]-1;
    roiExtDx = imcrop(extDx,roiBox);
    roiExtDy = imcrop(extDy,roiBox);
    % Now split the cell
    if isfield(prevStruct,'mesh') && size(prevStruct.mesh,1)>1
      mesh = prevStruct.mesh;
      if ismember(p.algorithm,[2 3])
        pcCell = splitted2model(mesh,p);
        pcCell = model2box(pcCell,roiBox,p.algorithm);
        pcCell = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,false,roiBox,thres,[frame cell]);
        pcCell = box2model(pcCell,roiBox,p.algorithm);
        cCell = model2geom(pcCell,p.algorithm);
      elseif ismember(p.algorithm,4)
        pcCell = align4IM(mesh,p);
        pcCell = model2box(pcCell,roiBox,p.algorithm);
        cCell = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,roiBox,thres,[frame cell]);
        cCell = box2model(cCell,roiBox,p.algorithm);
      end
      if isempty(cCell), gdisp(['Cell ' num2str(cell) ' refinement failed: error fitting shape']); continue; end
      mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth);
      if isempty(mesh), gdisp(['Cell ' num2str(cell) ' refinement failed: error getting mesh']); continue; end
      if size(pcCell,2)==1, model=cCell'; else model=cCell; end
      prevStruct.model = model;
      prevStruct.mesh = mesh;
      roiBox(1:2) = floor(max([min(min(mesh(:,[1 3]))) min(min(mesh(:,[2 4])))]-p.roiBorder,1));
      roiBox(3:4) = ceil(min([max(max(mesh(:,[1 3]))) max(max(mesh(:,[2 4])))]+p.roiBorder,[size(img,2) size(img,1)])-roiBox(1:2));
      prevStruct.box = roiBox;
      if ~isfield(prevStruct,'timelapse')
        if ismember(0,cellListN(1)==cellListN) || length(cellListN)<=1
          prevStruct.timelapse = 0;
        else
          prevStruct.timelapse = 1;
        end
      end
      cellList{frame}{cell} = getextradata(prevStruct);
      n = n+1;
    end
  end
  gdisp(['Refinement of ' num2str(n) ' cells succeeded'])
end

function lst = forcesplitcell(frame,cell,splitpos)
  % this function attempts to split the selected cell (it will add the
  % daughter cell to the cellList if successful, returns the numbers of
  % the two cells)
  
  global cellList cellListN p rawPhaseData
  
  lst = [];
  if checkparam(p,'invertimage','algorithm','erodeNum')
    gdisp('Splitting cells failed: one or more required parameters not provided.');
    return
  end
  if length(cellList)<frame || isempty(cellList{frame}) || length(cell)~=1 || isempty(cellList{frame}{cell}), return; end
  if frame>size(rawPhaseData,3), gdisp('Splitting cells failed: no image for this frame'); return; end
  if p.invertimage
    img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
  else
    img = rawPhaseData(:,:,frame);
  end
  imge = img2imge(img,p.erodeNum);%
  imge16 = img2imge16(img,p.erodeNum);
  thres = graythreshreg(imge,p.threshminlevel);
  bgr = phasebgr(imge,thres);
  [extDx,extDy] = getExtForces(imge,imge16,p);
  if isempty(extDx), gdisp('Force splitting cells failed: unable to get energy'); return; end
  prevStruct = cellList{frame}{cell};
  roiBox = prevStruct.box;%
  roiImg = imcrop(imge,roiBox);%
  roiExtDx = imcrop(extDx,roiBox);
  roiExtDy = imcrop(extDy,roiBox);
  
  % Now split the cell - copied from ProcessFrameI
  if ismember(p.algorithm,[2 3 4]) && isfield(prevStruct,'mesh') && size(prevStruct.mesh,1)>1
    mesh = prevStruct.mesh;
    if isempty(splitpos)
      res=isDivided(mesh-repmat(roiBox(1:2),size(mesh,1),2),roiImg,0,bgr);
      if isempty(res || res<0), lst = cell; gdisp('Cell division failed: no division at zero threshold'); return; end
    else
      res = splitpos;
    end
    mesh1 = flipud(mesh(res+1:end,:)); % daughter cell
    mesh2 = mesh(1:res-1,:); % mother cell
    if ismember(p.algorithm,[2 3])
      pcCell1 = splitted2model(mesh1,p);
      pcCell1 = model2box(pcCell1,roiBox,p.algorithm);
      pcCell1 = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell1,p,false,roiBox,thres,[frame cell]);
      pcCell1 = box2model(pcCell1,roiBox,p.algorithm);
      cCell1 = model2geom(pcCell1,p.algorithm);
      pcCell2 = splitted2model(mesh2,p);
      pcCell2 = model2box(pcCell2,roiBox,p.algorithm);
      pcCell2 = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell2,p,false,roiBox,thres,[frame cell]);
      pcCell2 = box2model(pcCell2,roiBox,p.algorithm);
      cCell2 = model2geom(pcCell2,p.algorithm);
    elseif p.algorithm==4
      pcCell1 = align4IM(mesh1,p);
      pcCell1 = model2box(pcCell1,roiBox,p.algorithm);
      cCell1 = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell1,p,roiBox,thres,[frame cell]);
      cCell1 = box2model(cCell1,roiBox,p.algorithm); % Corrected pcCell->cCell 2008/08/02
      pcCell2 = align4IM(mesh2,p);
      pcCell2 = model2box(pcCell2,roiBox,p.algorithm);
      cCell2 = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell2,p,roiBox,thres,[frame cell]);
      cCell2 = box2model(cCell2,roiBox,p.algorithm); % Corrected pcCell->cCell 2008/08/02
    end
    if isempty(cCell1) || isempty(cCell2), lst = cell; gdisp('Splitting cells failed: error fitting shape'); return; end

    % algorithm-specific portion of constructing the cell structures
    if size(pcCell1,2)==1, model=cCell1'; else model=cCell1; end
    cellStruct1.model = model;
    if size(pcCell2,2)==1, model=cCell2'; else model=cCell2; end
    cellStruct2.model = model;
    
  elseif p.algorithm==1
    getmesh = isfield(prevStruct,'mesh') && size(prevStruct.mesh,1)>1;
    if getmesh
      mesh = prevStruct.mesh;
      mask = poly2mask([mesh(:,1);flipud(mesh(:,3))]-roiBox(1),[mesh(:,2);flipud(mesh(:,4))]-roiBox(2),roiBox(4)+1,roiBox(3)+1);
    else
      contour = prevStruct.contour;
      mask = poly2mask(contour(:,1)-roiBox(1),contour(:,2)-roiBox(2),roiBox(4)+1,roiBox(3)+1);
    end
    [mask1,mask2]=splitonereg(roiImg.*mask);
    if isempty(mask1), gdisp(['Splitting cells failed: unable to split ' num2str(cell)]); return; end
    
    % Get and smooth the boundary of the region
    [ii,jj]=find(bwperim(mask1),1,'first');
    pp = bwtraceboundary(mask1,[ii,jj],'n',8,inf,'counterclockwise');
    fpp = frdescp(pp);
    cCell1 = ifdescp(fpp,p.fsmooth)+1; % +1?
    cCell1 = [cCell1;cCell1(1,:)]; % make the last point the same as the first one
    cCell1 = rot90(cCell1,2);
    cCell1(:,1) = cCell1(:,1)+roiBox(1)-1;
    cCell1(:,2) = cCell1(:,2)+roiBox(2)-1;
    [ii,jj]=find(bwperim(mask2),1,'first');
    pp=bwtraceboundary(mask2,[ii,jj],'n',8,inf,'counterclockwise');
    fpp = frdescp(pp);
    cCell2 = ifdescp(fpp,p.fsmooth)+1; 
    cCell2 = [cCell2;cCell2(1,:)]; % make the last point the same as the first one
    cCell2 = rot90(cCell2,2);
    cCell2(:,1) = cCell2(:,1)+roiBox(1)-1;
    cCell2(:,2) = cCell2(:,2)+roiBox(2)-1;
    
    % algorithm-specific posrtion of constructiong the cell structures
    cellStruct1.contour = cCell1;
    cellStruct2.contour = cCell2;
  end
  
  % finish constructiong the cell structures (same for all algorithms)
  cellStruct1.algorithm = p.algorithm;
  cellStruct2.algorithm = p.algorithm;
  if isfield(prevStruct,'mesh') && size(prevStruct.mesh,1)>1
    mesh1 = model2mesh(cCell1,p.meshStep,p.meshTolerance,p.meshWidth);
    mesh2 = model2mesh(cCell2,p.meshStep,p.meshTolerance,p.meshWidth);
    if isempty(mesh1) || isempty(mesh2), lst = cell; gdisp('Splitting cells failed: error getting mesh'); return; end
    cellStruct1.mesh = mesh1;
    cellStruct2.mesh = mesh2;
  end
  cellStruct1.stage = 1;
  cellStruct2.stage = 1;
  if ~isfield(prevStruct,'timelapse')
    if isempty(cellListN) || ismember(0,cellListN(1)==cellListN)
      prevStruct.timelapse = 0;
    else
      prevStruct.timelapse = 1;
    end
  end
  cellStruct1.timelapse = prevStruct.timelapse;
  cellStruct2.timelapse = prevStruct.timelapse;
  cellStruct1.box = roiBox;
  cellStruct2.box = roiBox;
  if ~p.forceindframes && prevStruct.timelapse
    if frame>1
      mdescendants = [];
      for cframe=frame-1:-1:1
        if length(cellList{cframe})>=cell && ~isempty(cellList{cframe}{cell})
          mdescendants = cellList{cframe}{cell}.descendants;
          break
        end
      end
      mdivisions = prevStruct.divisions(prevStruct.divisions~=frame);
      % The next block has been temporarily removed. The idea was 
      % that manual division on a later frame would assign the 
      % nexborn cell the same number as already produced cell. The 
      % way it is without the block, the new cell will get a new 
      % number, so after a division both progeny hes to be rerun on 
      % all the consequitive frames.
      if false%~isempty(mdescendants) && (mdescendants(end)<length(cellList{frame}) || isempty(cellList{frame}{mdescendants(end)}))
        daughter = mdescendants(end);
        cellStruct2.divisions = mdivisions;
        for cframe=1:frame-1
          if (length(cellList{cframe})>=daughter && ~isempty(cellList{cframe}{daughter})) ...
              || ( length(cellList{cframe})>=cell && ~isempty(cellList{cframe}{cell}) ...
                 && ismember(daughter,cellList{cframe}{cell}.descendants) )
            cellStruct1.divisions = cframe;
            cellStruct1.birthframe = cframe;
            break
          end
        end
      else
        birthframe = prevStruct.birthframe;
        mdivisions = [mdivisions frame];
        daughter = getdaughter(cell,sum(mdivisions~=birthframe),cellListN(frame));
        cellStruct2.divisions = mdivisions;
        cellStruct1.divisions = frame;
        cellStruct1.birthframe = frame;
      end
      cellStruct2.birthframe = prevStruct.birthframe;
      cellStruct1.ancestors = [prevStruct.ancestors(prevStruct.ancestors~=cell) cell];
      cellStruct2.ancestors = prevStruct.ancestors;
      cellStruct1.descendants = [];
      cellStruct2.descendants = [prevStruct.descendants(prevStruct.descendants~=daughter) daughter];
      cellStruct1.polarity = 1;
      cellStruct2.polarity = 1;
    else % case of frame==1
      daughter = length(cellList{1})+1;
      for lng = 1:length(cellList{1}) % try putting the daughter instead of a deleted cell to keep cellListN unchanged
        if isempty(cellList{1}{lng}), daughter=lng; break, end
      end
      if daughter>lng
        cellListN(1) = daughter;
        gdisp('The number of cells on the first frame (cellListN) changed: make sure to recompute all subsequent frames.')
      end
      cellStruct1.divisions = [];
      cellStruct2.divisions = [];
      cellStruct1.birthframe = 1;
      cellStruct2.birthframe = 1;
      cellStruct1.ancestors = [];
      cellStruct2.ancestors = [];
      cellStruct1.descendants = [];
      cellStruct2.descendants = [];
      cellStruct1.polarity = 0; % polarity is not set when the cells on the 1st frame are splitted
      cellStruct2.polarity = 0;
    end
  else
    % independent frames regime: keep the attributess of the mother
    % cell, set the attributes of the daughter cell as newly appeared
    daughter = length(cellList{frame})+1;
    cellListN(frame) = daughter;
    cellStruct1.divisions = [];
    cellStruct2.divisions = prevStruct.divisions;
    cellStruct1.ancestors = [];
    cellStruct2.ancestors = prevStruct.ancestors;
    cellStruct1.descendants = [];
    cellStruct2.descendants = prevStruct.descendants;
    cellStruct1.birthframe = frame;
    cellStruct2.birthframe = prevStruct.birthframe;
    cellStruct1.polarity = 0;
    cellStruct2.polarity = 0;
  end
  
  % save the data to cellList
  cellList{frame}{daughter} = cellStruct1;
  cellList{frame}{cell} = cellStruct2;
  cellList{frame}{daughter} = getextradata(cellList{frame}{daughter});
  cellList{frame}{cell} = getextradata(cellList{frame}{cell});
  lst = [cell daughter];
  if ~p.forceindframes && prevStruct.timelapse
    updatelineage(lst,frame)
  end
  
  % display the result message
  if ~p.forceindframes && prevStruct.timelapse
    gdisp(['Splitting cell ' num2str(cell) ' succeeded, saved as cells ' num2str(cell) ' and ' num2str(daughter) ', marked as independent.'])
  else
    gdisp(['Splitting cell ' num2str(cell) ' succeeded, saved as cells ' num2str(cell) ' and ' num2str(daughter) ', marked as mother & daughter.'])
  end
end

function updatelineage(lst,frame)
  % This function updates the lineage information for the group of
  % selected cells. It assumes that on the current frame is correct and
  % goes through the successive frames
  
  global cellList cellListN
  
  for i=1:length(lst)
    if length(cellList{frame})<lst(i) || isempty(cellList{frame}{lst(i)}) || ~isfield(cellList{frame}{lst(i)},'timelapse') ...
                      || ~cellList{frame}{lst(i)}.timelapse
      gdisp('Lineage could not be updated')
      return
    end
  end
  if length(cellListN)<=1 || ismember(0,cellListN(1)==cellListN) 
    gdisp('Lineage could not be updated')
    return
  end
  
  pcelllist = lst;
  for i=1:length(lst)
    pstructlist{i} = movestruct([],cellList{frame}{lst(1)});
  end
  for cframe=frame+1:length(cellList)
    ccelllist = [];
    cstructlist = {};
    for i=1:length(pcelllist)
      pcell = pcelllist(i);
      mstruct = pstructlist{i};
      daughter = getdaughter(pcell,length(mstruct.divisions),cellListN(cframe));
      if length(cellList{cframe})>=daughter && ~isempty(cellList{cframe}{daughter})
        mstruct.descendants = [mstruct.descendants daughter];
        mstruct.divisions = [mstruct.divisions cframe];
        dstruct.birthframe = cframe;
        dstruct.ancestors = [mstruct.ancestors pcell];
        dstruct.descendants = [];
        dstruct.divisions = cframe;
        ccelllist = [ccelllist pcell daughter];
        cstructlist = [cstructlist mstruct dstruct];
        cellList{cframe}{daughter} = movestruct(cellList{cframe}{daughter},dstruct);
      else
        ccelllist = [ccelllist pcell];
        cstructlist = [cstructlist mstruct];
      end
      if length(cellList{cframe})>=pcell && ~isempty(cellList{cframe}{pcell})
        cellList{cframe}{pcell} = movestruct(cellList{cframe}{pcell},mstruct);
      end
    end
    pcelllist = ccelllist;
    pstructlist = cstructlist;
  end
end
function s1=movestruct(s1,s2)
  % update lineage for the updatelineage function
  s1.birthframe = s2.birthframe;
  s1.ancestors = s2.ancestors;
  s1.descendants = s2.descendants;
  s1.divisions = s2.divisions;
end

function [mask1,mask2] = splitonereg(img)
  % this function splits a region along the inverse drainage divide 
  % containing the lowest pass between any two peaks inside of a region
  
  % find watersheds of the inverted image
  img = im2uint16(img);
  mask = img~=0;
  wshed = mask.*double(watershed(2^16-1-img));
  
  % find the highest level, which leaves all wsheds
  clevel = 2^16;
  level1 = 0;
  nwsheds = max(max(wshed));
  for i=1:16
    clevel = clevel/2;
    cimg = (img>clevel+level1).*wshed;
    if max(max(bwlabel(cimg>0)))==nwsheds
      level1 = level1+clevel;
    end
  end
  
  % find the highest level below previous, which keeps one region
  clevel = 2^16;
  level2 = 0;
  for i=1:16
    clevel = clevel/2;
    if level2+clevel<level1
      cimg = (img>clevel+level2);
      if max(max(bwlabel(cimg>0)))==1
        level2 = level2+clevel;
      end
    end
  end
  
  % get the masks
  cimg = img>level2+1;
  lbl = bwlabel(cimg|(wshed>0),4);
  [row,col] = find(mask.*(lbl==0));
  newmask = lbl;
  for i=1:length(row)
    ind1 = max(row(i)-1,1):min(row(i)+1,size(lbl,1));
    ind2 = max(col(i)-1,1):min(col(i)+1,size(lbl,2));
    img2 = reshape(img(ind1,ind2),1,[]);
    lbl2 = reshape(lbl(ind1,ind2),1,[]);
    [v,ind]=max(img2);
    if sum(lbl2(ind)==1)>0 && sum(lbl2(ind)>1)==0
      newmask(row(i),col(i))=1;
    elseif sum(lbl2(ind)>1)>0 && sum(lbl2(ind)==1)==0
      newmask(row(i),col(i))=2;
    end
  end
  mask1 = newmask==1;
  mask2 = newmask>1; % note: mask2 may contain more than one region!
  if max(max(mask2))==0, mask1=[]; mask2=[]; end
end

function res = forcejoincells(frame,lst) % TODO: write for algorithm 1
  % this function tries to join all the cells in a given group
  % works for multiple cells, has an extra parameter for dilation
  % currently only works algorithms 2-4
  
  global cellList cellListN p se rawPhaseData maskdx maskdy
  res = lst; % output the input list is failed
  if checkparam(p,'invertimage','algorithm','erodeNum','roiBorder','meshStep','meshTolerance','meshWidth','joindilate')
    gdisp('Joining cells failed: one or more required parameters not provided.');
    return
  end  
  if ~ismember(p.algorithm,[2 3 4]), gdisp('Joining cells failed: joining routine is only available for algorithms 2-4.'); return; end
  if length(cellList)<frame || isempty(cellList{frame}) || length(lst)<1, return; end
  if isempty(cellList{frame}{lst(1)}) || isempty(cellList{frame}{lst(2)}), return; end
  pstruct = cellList{frame}{min(lst)};
  if ~isfield(pstruct,'timelapse')
    if length(cellListN)<=1 || ismember(0,cellListN(1)==cellListN)
      pstruct.timelapse = 0;
    else
      pstruct.timelapse = 1;
    end
  end
  if ~p.forceindframes && pstruct.timelapse && length(lst)~=2
    gdisp('Joining cells failed: in the timelapse mode, exactly two cells must be selected.')
    return; 
  end
  
  if frame>size(rawPhaseData,3), return; end
  if p.invertimage
    img = max(max(max(rawPhaseData)))-rawPhaseData(:,:,frame);
  else
    img = rawPhaseData(:,:,frame);
  end
  imge = img2imge(img,p.erodeNum);
  imge16 = img2imge16(img,p.erodeNum);
  thres = graythreshreg(imge,p.threshminlevel);
  [extDx,extDy] = getExtForces(imge,imge16,p);
  if isempty(extDx), gdisp('Force joining cells failed: unable to get energy'); return; end
  
  % Start the joining procedire immediately without checking distance
  border = p.roiBorder;
  mesh = {};
  box = zeros(length(lst),4);
  for i=1:length(lst)
    mesh{i} = cellList{frame}{lst(i)}.mesh;
    box(i,:) = [min(min(mesh{i}(:,[1 3]))) min(min(mesh{i}(:,[2 4]))) (-max(max(mesh{i}(:,[1 3])))) (-max(max(mesh{i}(:,[2 4]))))];
  end
  roiBox = min(box);
  roiBox = [max(floor(roiBox(1:2)-border),1) min(ceil(-roiBox(3:4)+border),[size(img,2) size(img,1)])];
  roiBox = [roiBox(1:2) roiBox(3:4)-roiBox(1:2)]+1;
  
  % Create a mask
  mask = poly2mask([mesh{1}(:,1);flipud(mesh{1}(:,3))]-roiBox(1)+1,[mesh{1}(:,2);flipud(mesh{1}(:,4))]-roiBox(2)+1,roiBox(4),roiBox(3));
  if p.joindilate<0
    for j=1:-p.joindilate
      mask = imerode(mask,se);
    end
  end
  p1 = mesh{1}(1,:);
  p2 = mesh{1}(end,:);
  p1a = mesh{1}(4,:);
  p2a = mesh{1}(end-3,:);
  lst2 = 2:length(lst);
  for i=1:length(lst)-1
    dsn = [];
    ind = [];
    for j=lst2
      [dsn2,ind2] = min([dstm(p1-mesh{j}(1,:)) dstm(p1-mesh{j}(end,:)) dstm(p2-mesh{j}(1,:)) dstm(p2-mesh{j}(end,:))]);
      dsn = [dsn dsn2];
      ind = [ind ind2];
    end
    [dsn3,ind3] = min(dsn);
    ind4 = ind(ind3);
    mesh2 = mesh{lst2(ind3)};
    lst2 = lst2(lst2~=lst2(ind3));
    if ind4==1 || ind4==2, p1b = p1a; else p1b = p2a; end
    if ind4==1 || ind4==3, p2b = mesh2(4,:); else p2b = mesh2(end-3,:); end
    mask2 = poly2mask([mesh2(:,1);flipud(mesh2(:,3))]-roiBox(1)+1,[mesh2(:,2);flipud(mesh2(:,4))]-roiBox(2)+1,roiBox(4),roiBox(3));
    mask3 = poly2mask([p1b(1) p2b(1) p2b(3) p1b(3)]-roiBox(1)+1,[p1b(2) p2b(2) p2b(4) p1b(4)]-roiBox(2)+1,roiBox(4),roiBox(3));
    mask4 = poly2mask([p1b(1) p2b(3) p2b(1) p1b(3)]-roiBox(1)+1,[p1b(2) p2b(4) p2b(2) p1b(4)]-roiBox(2)+1,roiBox(4),roiBox(3));
    if p.joindilate<0
      for j=1:-p.joindilate
        mask2 = imerode(mask2,se);
        mask3 = imerode(mask3,se);
        mask4 = imerode(mask4,se);
      end
    end
    if sum(sum(mask.*mask2))>100, mask3=0; mask4=0; end
    mask = min(1,mask+mask2+mask3+mask4);
    if ind4==1, p1 = mesh2(end,:); p1a = mesh2(end-3,:); end
    if ind4==2, p1 = mesh2(1,:); p1a = mesh2(4,:); end
    if ind4==3, p2 = mesh2(end,:); p2a = mesh2(end-3,:); end
    if ind4==4, p2 = mesh2(1,:); p2a = mesh2(4,:); end
  end
  if p.joindilate>0
    for j=1:p.joindilate
      mask = imdilate(mask,se);
    end
  end
  edg = bwperim(mask);

  pmap = 1 - edg;
  f1=true;
  while f1
    pmap1 = 1 - edg + imerode(pmap,se);
    f1 = max(max(pmap1-pmap))>0;
    pmap = pmap1;
  end;
  pmapEnergy = pmap + 0.1*pmap.^2;
  pmapDx = imfilter(pmapEnergy,maskdx); % distance forces
  pmapDy = imfilter(pmapEnergy,maskdy); 
  pmapDxyMax = 10;
  pmapDx = pmapDx/pmapDxyMax; % normalize to make the max force equal to 1
  pmapDy = pmapDy/pmapDxyMax;

  % roiBox2 = roiBox+[1 1 0 0];%[roiBox(2) roiBox(1) roiBox(4)-1 roiBox(3)-1]; % standard format
  roiImg = imcrop(imge,roiBox);
  roiExtDx = imcrop(extDx,roiBox);
  roiExtDy = imcrop(extDy,roiBox);
  roiBox(3:4) = [size(roiExtDx,2) size(roiExtDx,1)]-1; %!

  prop = regionprops(bwlabel(mask),'orientation','centroid'); % TODO: Check
  theta = prop(1).Orientation*pi/180;
  x0 = prop(1).Centroid(1);
  y0 = prop(1).Centroid(2);

  if p.algorithm==2
    pcCell0 = [theta;x0;y0;zeros(p.Nkeep+1,1)];
  elseif p.algorithm==3
    pcCell0 = [x0;y0;theta;0;zeros(p.Nkeep+1,1)];
  end
  if ismember(p.algorithm,[2 3])
    [pcCell,fitquality] = align(mask,pmapDx,pmapDy,pmapDx*0,pcCell0,p,true,roiBox,0.5,[frame lst(1)]);
    [pcCell,fitquality] = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,false,roiBox,thres,[frame lst(1)]);
  elseif p.algorithm == 4
    pcCell = align4I(mask,p);
    [pcCell,fitquality] = align4(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,roiBox,thres,[frame lst(1)]);
  end
  pcCell = box2model(pcCell,roiBox,p.algorithm);
  cCell = model2geom(pcCell,p.algorithm);
  mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth);

  if length(mesh)<5
    gdisp(['Tried to join cells ' num2str(lst) ' - failed, resplit again.'])
  else
    % construct the lineage data
    if ~p.forceindframes && pstruct.timelapse && frame>1
      mcell = min(lst);
      dcell = lst(lst~=mcell);
      if frame>1
        dcellexist = false;
        for cframe=frame-1:-1:1
          if mcell<=length(cellList{cframe}) && ~isempty(cellList{cframe}{mcell})
            if ismember(dcell,cellList{cframe}{mcell}.descendants)
              dcellexist = true;
              break
            end
          end
        end
        if ~dcellexist
          f = 1+length(pstruct.divisions)-find(pstruct.descendants==dcell);
          pstruct.descendants = pstruct.descendants(pstruct.descendants~=dcell);
          pstruct.divisions(f) = [];
          cellList{frame}{mcell}.divisions = pstruct.divisions;
          cellList{frame}{mcell}.descendants = pstruct.descendants;
          cellList{frame}{mcell}.timelapse = 1;
          cellList{frame}{dcell} = [];
          updatelineage(mcell,frame);
        end
      end
    elseif ~p.forceindframes && pstruct.timelapse && frame==1
      mcell = min(lst);
      dcell = lst(lst~=mcell);
      cellList{frame}{mcell}.timelapse = 1;
    else
      mcell = lst(1);
      dcell = lst(2:end);
    end
    
    % save the data to cellList
    cellList{frame}{mcell}.mesh = mesh;
    if size(pcCell,2)==1, model=reshape(pcCell,[],1); else model=pcCell; end
    cellList{frame}{mcell}.model = model; % TODO: check if works for alg. 2-3
    cellList{frame}{mcell}.box = roiBox;
    for i=dcell, cellList{frame}{i} = []; end
    cellList{frame}{mcell}=getextradata(cellList{frame}{mcell});
    gdisp(['Joined cells ' num2str(lst) ' - success, saved as cell ' num2str(mcell) '.'])
    res = mcell; % output the number of the new cell is succeeded
  end
  
  function res = dstm(a)
    res = (a(1)+a(3))^2+(a(2)+a(4))^2;
  end
end

function [pcCell,cCell] = splitted2model(mesh,p)
% This function fits the model to a half of a splitted cell defined by mesh
% such that the Nth model point corresponds to the top row of mesh
if size(mesh,1)<4, pcCell = -1; cCell = -1; return; end
mesh = [repmat(([1;2;3]*mesh(1,1:2)+[3;2;1]*mesh(1,3:4))/4,1,2);mesh];
global se maskdx maskdy;
border = p.roiBorder;
roiBox = round([min(min(mesh(:,[1 3])))-border min(min(mesh(:,[2 4])))-border...
        max(max(mesh(:,[1 3])))-min(min(mesh(:,[1 3])))+2*border...
        max(max(mesh(:,[2 4])))-min(min(mesh(:,[2 4])))+2*border]);
edg = zeros(roiBox(3:4)+1);
edg(sub2ind(roiBox(3:4)+1,round([mesh(:,1);mesh(:,3)]-roiBox(1)),round([mesh(:,2);mesh(:,4)]-roiBox(2)))) = 1;
edg = edg';

pmap = 1 - edg;
f1=true;
while f1
  pmap1 = 1 - edg + imerode(pmap,se);
  f1 = max(max(pmap1-pmap))>0;
  pmap = pmap1;
end;

pmapEnergy = pmap + 0.1*pmap.^2;
pmapDx = imfilter(pmapEnergy,maskdx); % distance forces
pmapDy = imfilter(pmapEnergy,maskdy); 
pmapDxyMax = 10;
pmapDx = pmapDx/pmapDxyMax; % normalize to make the max force equal to 1
pmapDy = pmapDy/pmapDxyMax;

mask = poly2mask([mesh(:,1);flipud(mesh(:,3))]-roiBox(1),[mesh(:,2);flipud(mesh(:,4))]-roiBox(2),roiBox(4)+1,roiBox(3)+1);

for i=1%:2
  theta = pi-angle(mesh(end,1)+mesh(end,3)-mesh(1,1)-mesh(1,3)+...
        j*(mesh(end,2)+mesh(end,4)-mesh(1,2)-mesh(1,4)));
  x0 = mean(mean(mesh(:,[1 3])))-roiBox(1);
  y0 = mean(mean(mesh(:,[2 4])))-roiBox(2);
  
  if p.algorithm==2
    pcCell0 = [theta;x0;y0;zeros(p.Nkeep+1,1)];
  elseif p.algorithm==3
    pcCell0 = [x0;y0;theta;0;zeros(p.Nkeep+1,1)];
  elseif p.algorithm==4
    %pcCell0 = align4I(mask,p);
  end

  % Making first variant of the model
  if ismember(p.algorithm,[2 3])
    % first approximation - to the exact shape of the selected region
    [pcCell,fitquality] = align(mask,pmapDx,pmapDy,pmapDx*0,pcCell0,p,true,roiBox,0.5,[0 0 1]);
    % adjustment of the model to the external energy map
    %[pcCell,fitquality] = align(roiImg,roiExtDx,roiExtDy,roiExtDx*0,pcCell,p,false,roiBox,thres,[0 0]);
  elseif p.algorithm == 4
    pcCell = align4I(mask,p);
    fitquality = 0;
    %[pcCell,fitquality] = align4(roiImg,roiExtDx,roiExtDy,pmap*0,pcCell,p,roiBox,thres,[0 0]);
  end
  %gdisp(['fitquality aligning = ' num2str(fitquality)])
      
  % converting from box to global coordinates
  pcCell = box2model(pcCell,roiBox,p.algorithm);
  % obtaining the shape of the cell in geometrical representation
  cCell = model2geom(pcCell,p.algorithm);
  if isempty(cCell), return; end

  if (cCell(end,1)-mean(mesh(1,[1 3])))^2+(cCell(end,2)-mean(mesh(1,[2 4])))^2 >...
    (cCell(end,1)-mean(mesh(end,[1 3])))^2+(cCell(end,2)-mean(mesh(end,[2 4])))^2
    if p.algorithm==2
      pcCell(1) = pcCell(1) + pi;
    elseif p.algorithm==3
      pcCell(3) = pcCell(3) + pi;
    end
  else
    break;
  end
end
gdisp(['splitting fitquality = ' num2str(fitquality)])
end



% ------------- Functions, that used to be in separate files --------------

function res = model2mesh(coord,stp,tolerance,lng)
% This function performs a medial axis transform to a non-branching axis
% Takes the outline coordinates, step size on the final centerline,
% tolerance to non-centrality, and the length of the ribs. The last
% parameter should be longer than the ribs, but should not be too long to
% intersect the countour agan: though most cases of >1 intersection will
% be resolved by the function, some will not. Outputs the coordinates of
% the centerline points.
  delta = 1E-10;
  
  % voronoi transform
  while true
    if length(coord)<=2, res=0; return; end
    if abs(coord(1,1)-coord(end,1))<0.00001 && abs(coord(1,2)-coord(end,2))<0.00001
      coord = coord(1:end-1,:);
    else
      break
    end
  end
  warning('off','MATLAB:delaunay:DuplicateDataPoints')
  warning('off','MATLAB:TriRep:PtsNotInTriWarnId')
  [vx,vy] = voronoi(coord(:,1),coord(:,2));
  warning('on','MATLAB:delaunay:DuplicateDataPoints')
  warning('on','MATLAB:TriRep:PtsNotInTriWarnId')
  
  % remove vertices crossing the boundary
  vx = reshape(vx,[],1);
  vy = reshape(vy,[],1);
  q = intxy2(vx,vy,coord([1:end 1],1),coord([1:end 1],2));
  vx = reshape(vx(~[q;q]),2,[]);
  vy = reshape(vy(~[q;q]),2,[]);
  
  % remove vertices outside
  q = logical(inpolygon(vx(1,:),vy(1,:),coord(:,1),coord(:,2))...
       .* inpolygon(vx(2,:),vy(2,:),coord(:,1),coord(:,2)));
  vx = reshape(vx([q;q]),2,[]);
  vy = reshape(vy([q;q]),2,[]);
  
  % remove isolated points
  if isempty(vx), res=0;return; end
  t = ~((abs(vx(1,:)-vx(2,:))<delta)&(abs(vy(1,:)-vy(2,:))<delta));
  vx = reshape(vx([t;t]),2,[]);
  vy = reshape(vy([t;t]),2,[]);
     
  % remove branches
  vx2=[];
  vy2=[];
  while true
    for i=1:size(vx,2)
      if((sum(sum((abs(vx-vx(1,i))<delta)&(abs(vy-vy(1,i))<delta)))>1)&&(sum(sum((abs(vx-vx(2,i))<delta)&(abs(vy-vy(2,i))<delta)))>1))
        vx2=[vx2 vx(:,i)];
        vy2=[vy2 vy(:,i)];
      end
    end
    if size(vx,2)-size(vx2,2)<=2,
      vx3 = vx2;
      vy3 = vy2;
      break;
    else
      vx = vx2;
      vy = vy2;
      vx2 = [];
      vy2 = [];
    end
  end
  vx3 = vx;
  vy3 = vy;
  
  % check that there are no cycles
  vx2 = [];
  vy2 = [];
  while size(vx,2)>1
    for i=1:size(vx,2)
      if((sum(sum((abs(vx-vx(1,i))<delta)&(abs(vy-vy(1,i))<delta)))>1)&&(sum(sum((abs(vx-vx(2,i))<delta)&(abs(vy-vy(2,i))<delta)))>1))
        vx2=[vx2 vx(:,i)];
        vy2=[vy2 vy(:,i)];
      end
    end
    if size(vx,2)-size(vx2,2)<=1
      res = 0;
      return;
    else
      vx = vx2;
      vy = vy2;
      vx2 = [];
      vy2 = [];
    end
  end
  vx = vx3;
  vy = vy3;
  if isempty(vx) || size(vx,1)<2, res=0;return;end

  % sort points
  vx2=[];
  vy2=[];
  for i=1:size(vx,2) % in this cycle find the first point
    if sum(sum(abs(vx-vx(1,i))<delta & abs(vy-vy(1,i))<delta))==1
      vx2=vx(:,i)';
      vy2=vy(:,i)';
      break
    elseif sum(sum(abs(vx-vx(2,i))<delta & abs(vy-vy(2,i))<delta))==1
      vx2=fliplr(vx(:,i)');
      vy2=fliplr(vy(:,i)');
      break
    end
  end
  k=2;
  while true % in this cycle sort all points after the first one
    f1=find(abs(vx(1,:)-vx2(k))<delta & abs(vy(1,:)-vy2(k))<delta & (abs(vx(2,:)-vx2(k-1))>=delta | abs(vy(2,:)-vy2(k-1))>=delta));
    f2=find(abs(vx(2,:)-vx2(k))<delta & abs(vy(2,:)-vy2(k))<delta & (abs(vx(1,:)-vx2(k-1))>=delta | abs(vy(1,:)-vy2(k-1))>=delta));
    if f1>0
      vx2 = [vx2 vx(2,f1)];
      vy2 = [vy2 vy(2,f1)];
    elseif f2>0
      vx2 = [vx2 vx(1,f2)];
      vy2 = [vy2 vy(1,f2)];
    else
      break
    end
    k=k+1;
  end

  skel=[vx2' vy2'];
  
  if size(vx2,2)<=1
    res = 0;
    return;
  end

  % interpolate skeleton to equal step, extend outside of the cell and smooth
  % tolerance=0.001;
  d=diff(skel,1,1);
  l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
  if l(end)>=stp
    skel = [interp1(l,vx2,0:stp:l(end))' interp1(l,vy2,0:stp:l(end))'];
  else
    skel = [vx2(1) vy2(1);vx2(end) vy2(end)];
  end
  if size(skel,1)<=1 % || length(skel)<4
    res = 0;
    return;
  end
  lng0 = l(end);
  sz = lng0/stp;
  L = size(coord,1);
  coord = [coord;coord(1,:)];
  % lng = 100;
  skel2 = [skel(1,:)*(lng/stp+1) - skel(2,:)*lng/stp; skel;...
       skel(end,:)*(lng/stp+1) - skel(end-1,:)*lng/stp];
  d=diff(skel2,1,1);
  l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
  [l,i] = unique(l);
  skel2 = skel2(i,:);
  if length(l)<2 || size(skel2,1)<2, res=0; return; end
  % find the intersection of the 1st skel2 segment with the contour, the
  % closest to one of the poles (which will be called 'prevpoint')
  [tmp1,tmp2,indS,indC]=intxyMulti(skel2(2:-1:1,1),skel2(2:-1:1,2),coord(:,1),coord(:,2));
  [tmp1,prevpoint] = min([min(modL(indC,1)) min(modL(indC,L/2+1))]);
  if prevpoint==2, prevpoint=L/2; end
  % prevpoint = mod(round(indC(1))-1,L)+1;
  skel3=spsmooth(l,skel2',tolerance,0:stp:l(end))';%1:stp:

  % recenter and smooth again the skeleton
  [pintx,pinty,q]=skel2mesh(skel3);
  if length(pintx)<sz-1
    skel3=spsmooth(l,skel2',tolerance/100,0:stp:l(end))';
    [pintx,pinty,q]=skel2mesh(skel3);
  end
  if ~q || length(pintx)<sz-1  || length(skel3)<4
    res=0;
    return
  end
  skel = [mean(pintx,2) mean(pinty,2)];
  d=diff(skel,1,1);
  l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
  [l,i] = unique(l);
  skel = skel(i,:);
  if length(l)<2 || size(skel,1)<2, res=0; return; end
  skel=spsmooth(l,skel',tolerance,-3*stp:stp:l(end)+4*stp)';
  
  % get the mesh
  [pintx,pinty,q]=skel2mesh(skel);
  if ~q 
    res=0;
    return
  end
  res = [pintx(:,1) pinty(:,1) pintx(:,2) pinty(:,2)];
  if (pintx(1,1)-coord(end,1))^2+(pinty(1,1)-coord(end,2))^2 > (pintx(end,1)-coord(end,1))^2+(pinty(end,1)-coord(end,2))^2
    res = flipud(res);
  end
  if numel(res)==1 || length(res)<=4, res=0; gdisp('Unable to create mesh'); end
  if length(res)>1 && (res(1,1)~=res(1,3) || res(end,1)~=res(end,3)), res=0; gdisp('Mesh creation error! Cell rejected'); end

function out = modL(in,shift)
  out = mod(in-shift,L);
  out = min(out,L-out);
end
  
function [pintx,pinty,q]=skel2mesh(sk)
  % This function finds intersections of ribs with the contour
  % To be used in "model2mesh" function
  
  if isempty(sk), pintx=[]; pinty=[]; q=false; return; end
  % Find the intersection of the skel with the contour closest to prevpoint
  pintx=[];
  pinty=[];
  [intX,intY,indS,indC]=intxyMulti(sk(:,1),sk(:,2),coord(:,1),coord(:,2));
  if isempty(intX) || isempty(indC) || isempty(prevpoint)
    q=false;
    return;
  end
  [prevpoint,ind] = min(modL(indC,prevpoint));
  prevpoint = indC(ind);
  indS=indS(ind);
  if indS>(size(sk,1)+1-indS)
    sk = sk(ceil(indS):-1:1,:);
  else
    sk = sk(floor(indS):end,:);
  end
  % 2. define the first pair of intersections as this point
  % 3. get the list of intersections for the next pair
  % 4. if more than one, take the next in the correct direction
  % 5. if no intersections found in the reqion between points, stop
  % 6. goto 3.
  % Define the lines used to compute intersections
  d=diff(sk,1,1);
  plinesx1 = repmat(sk(1:end-1,1),1,2)+lng/stp*d(:,2)*[0 1];
  plinesy1 = repmat(sk(1:end-1,2),1,2)-lng/stp*d(:,1)*[0 1];
  plinesx2 = repmat(sk(1:end-1,1),1,2)+lng/stp*d(:,2)*[0 -1];
  plinesy2 = repmat(sk(1:end-1,2),1,2)-lng/stp*d(:,1)*[0 -1];
  % Allocate memory for the intersection points
  pintx = zeros(size(sk,1)-1,2);
  pinty = zeros(size(sk,1)-1,2);
  % Define the first pair of intersections as the prevpoint
  pintx(1,:) = [intX(ind) intX(ind)];
  pinty(1,:) = [intY(ind) intY(ind)];
  prevpoint1 = prevpoint;
  prevpoint2 = prevpoint;
  
  % for i=1:size(d,1), plot(plinesx(i,:),plinesy(i,:),'r'); end % Testing
  q=true;
  fg = 1;
  jmax = size(sk,1)-1;
  for j=2:jmax
    % gdisp(['Use 1: ' num2str(size(plinesx1(j,:))) ' ' num2str(size(coord(:,1))) ' j=' num2str(j)]); % Testing
    [pintx1,pinty1,tmp,indC1]=intxyMulti(plinesx1(j,:),plinesy1(j,:),coord(:,1),coord(:,2),floor(prevpoint1),1);%
    [pintx2,pinty2,tmp,indC2]=intxyMulti(plinesx2(j,:),plinesy2(j,:),coord(:,1),coord(:,2),ceil(prevpoint2),-1);%
    if (~isempty(pintx1))&&(~isempty(pintx2))
      if pintx1~=pintx2
        if fg==3
          break;
        end
        fg = 2;
        [prevpoint1,ind1] = min(modL(indC1,prevpoint1));
        [prevpoint2,ind2] = min(modL(indC2,prevpoint2));
        prevpoint1 = indC1(ind1); 
        prevpoint2 = indC2(ind2);
        pintx(j,:)=[pintx1(ind1) pintx2(ind2)];
        pinty(j,:)=[pinty1(ind1) pinty2(ind2)];
      else
        q=false;
        return
      end
    elseif fg==2
      fg = 3;
    end
  end
  pinty = pinty(pintx(:,1)~=0,:);
  pintx = pintx(pintx(:,1)~=0,:);
  [intX,intY,indS,indC]=intxyMulti(sk(:,1),sk(:,2),coord(:,1),coord(:,2));
  [prevpoint,ind] = max(modL(indC,prevpoint));
  pintx = [pintx;[intX(ind) intX(ind)]];
  pinty = [pinty;[intY(ind) intY(ind)]];
  nonan = ~isnan(pintx(:,1))&~isnan(pinty(:,1))&~isnan(pintx(:,2))&~isnan(pinty(:,2));
  pintx = pintx(nonan,:);
  pinty = pinty(nonan,:);
end
end

function res=intxy2(ax,ay,bx,by)
% modified intxy, finds the first point along line a where an intersection 
% rewritten to accept unsorted sets, 2xN of N segments
% finds all intersections, only gives the row numbers in the first set
  if isempty(ax)
    res=[];
  else
    res=intxy2C(ax,ay,bx,by);
  end
% below is a MATLAB version of the C/C++ routine above, slow but same results
% 
%   res = zeros(size(ax,2),1);
%   for i=1:size(ax,2)
%     %vector to next vertex in a
%     u=[ax(2,i)-ax(1,i) ay(2,i)-ay(1,i)];
%     %go through each vertex in b
%     for j=1:size(bx,2)
%       %check for intersections at the vortices
%       if (ax(1,i)==bx(1,j) && ay(1,i)==by(1,j)) || (ax(2,i)==bx(1,j) && ay(2,i)==by(1,j))...
%         || (ax(1,i)==bx(2,j) && ay(1,i)==by(2,j)) ||(ax(2,i)==bx(2,j) && ay(2,i)==by(2,j))
%         res(i) = 1;
%         continue
%       end
%       %vector from ith vertex in a to jth vertex in b
%       v=[bx(2,j)-ax(1,i) by(2,j)-ay(1,i)];
%       %vector from ith+1 vertex in a to jth vertex in b
%       w=[bx(1,j)-ax(2,i) by(1,j)-ay(2,i)];
%       %vector from ith vertex of a to jth-1 vertex of b
%       vv=[bx(1,j)-ax(1,i) by(1,j)-ay(1,i)];
%       %vector from jth-1 vertex of b to jth vertex of b
%       z=[bx(2,j)-bx(1,j) by(2,j)-by(1,j)];
%       %cross product of u and v
%       cpuv=u(1)*v(2)-u(2)*v(1);
%       %cross product of u and vv
%       cpuvv=u(1)*vv(2)-u(2)*vv(1);
%       %cross product of v and z
%       cpvz=v(1)*z(2)-v(2)*z(1);
%       %cross product of w and z
%       cpwz=w(1)*z(2)-w(2)*z(1);
%       % check if there is an intersection
%       if cpuv*cpuvv<0 && cpvz*cpwz<0
%         res(i) = 1;
%       end
%     end
%   end
end

function res = graythreshreg(img,varargin)
    % threshold calculated in a regionSelectionRect region
    if ~isempty(varargin), flevel=varargin{1}; else flevel=0; end
    global regionSelectionRect
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
            %b = graythresh(c(c>=level));
            %EDIT SH 130204
            b = myMultiOtsu(c(c>=level));
        else
            %b = graythresh(a);
            %EDIT SH 130204 
            b = myMultiOtsu(a);
        end
        function b = myMultiOtsu(a)
           %EDIT SH 130204 - 2 level otsu to allow for phase contrast haloing
           [dum1, dum2 , thres] = multiOtsu(a,3);
           %multi otsu assumes 8 bit im. mut convert to imdouble range
           thres = double(thres)/255;
           %want the brightest level
           b = thres(2);
        end
    end
end

%%
%% Shifting frames

function [shiftX,shiftY]=alignframes(A,depth)
  
  mrg = round(min(size(A,1),size(A,2))*0.05);
  fld = [size(A,1)-2*mrg size(A,2)-2*mrg];
  time1 = clock;
  x0=[0 1 1 0 -1 -1 -1 0 1];
  y0=[0 0 1 1 1 0 -1 -1 -1];
  nframes = size(A,3);
  %field = [size(A,1) size(A,2)];
  memory = double(A(:,:,1));
  score = zeros(1,9);
  shiftX = zeros(1,nframes);
  shiftY = zeros(1,nframes);
  for frame=2:nframes
    memory = memory*(1-1/depth) + double(A(:,:,frame-1))/depth;
    sc2memory = imresize(memory,0.5);
    sc4memory = imresize(memory,0.25);
    cframe = double(A(:,:,frame));
    sc2cframe = imresize(cframe,0.5);
    sc4cframe = imresize(cframe,0.25);
    [xF,yF] = alignoneframe(sc4cframe,sc4memory,0,0,round(mrg/4));
    [xF,yF] = alignoneframe(sc2cframe,sc2memory,2*xF,2*yF,round(mrg/2));
    [xF,yF] = alignoneframe(cframe,memory,2*xF,2*yF,mrg);
    shiftX(frame) = xF;
    shiftY(frame) = yF;
    if depth>1 && (xF~=0 || yF~=0)
      fldX = max(1,xF+1):min(fld(1),fld(1)+xF);
      fldY = max(1,yF+1):min(fld(2),fld(2)+yF);
      fldX2 = max(1,-xF+1):min(fld(1),fld(1)-xF);
      fldY2 = max(1,-yF+1):min(fld(2),fld(2)-yF);
      memory(fldX2,fldY2) = memory(fldX,fldY);
    end
    gdisp(['frame = ' num2str(frame) ', shift ' num2str(xF) ',' num2str(yF) ' pixels'])
  end
  shiftX = cumsum(shiftX);
  shiftY = cumsum(shiftY);
  time2 = clock;
  gdisp(['Finised, elapsed time ' num2str(etime(time2,time1)) ' s']);  

  function [x,y] = alignoneframe(cframe,memory,x,y,margin)
    cframetmp = cframe(margin+1:end-margin,margin+1:end-margin);
    memorytmp = memory(margin+1:end-margin,margin+1:end-margin);
    field = [size(cframe,1)-2*margin size(cframe,2)-2*margin];
    while true
      for j=1:9
        dx=x0(j);
        dy=y0(j);
        xJ = x + dx;
        yJ = y + dy;
        fieldX = max(1,xJ+1):min(field(1),field(1)+xJ);
        fieldY = max(1,yJ+1):min(field(2),field(2)+yJ);
        fieldX2 = max(1,-xJ+1):min(field(1),field(1)-xJ);
        fieldY2 = max(1,-yJ+1):min(field(2),field(2)-yJ);
        score(j) = corel(memorytmp(fieldX,fieldY),cframetmp(fieldX2,fieldY2));
      end
      [scmax,ind] = max(score);
      if ind==1, break; end
      x = x+x0(ind);
      y = y+y0(ind);
    end
  end
end

function B = shiftstack(A,varargin)
  if length(varargin)==1
    shift=varargin{1};
  else
    shift.x=varargin{1};
    shift.y=varargin{2};
  end
  if ~iscell(A)
    field = [size(A,1) size(A,2)];
    B = ones(size(A),class(A));
    for frame = 1:size(A,3)
      xJ = shift.x(frame);
      yJ = shift.y(frame);
      fieldX = max(1,xJ+1):min(field(1),field(1)+xJ);
      fieldY = max(1,yJ+1):min(field(2),field(2)+yJ);
      fieldX2 = max(1,-xJ+1):min(field(1),field(1)-xJ);
      fieldY2 = max(1,-yJ+1):min(field(2),field(2)-yJ);
      B(:,:,frame) = B(:,:,frame)*max(max(A(:,:,frame)));
      B(fieldX,fieldY,frame) = A(fieldX2,fieldY2,frame);
    end
  else
    B = A;
    for frame = 1:min(length(A),length(shift.x))
      xJ = shift.x(frame);
      yJ = shift.y(frame);
      if length(A)>=frame && ~isempty(A{frame})
        for cell=1:length(A{frame})
          if ~isempty(A{frame}{cell}) && isfield(A{frame}{cell},'mesh') && length(A{frame}{cell}.mesh)>1
            B{frame}{cell}.mesh(:,[1 3]) = A{frame}{cell}.mesh(:,[1 3])+yJ;
            B{frame}{cell}.mesh(:,[2 4]) = A{frame}{cell}.mesh(:,[2 4])+xJ;
          end
        end
      end
    end
  end
end

function y=corel(X,Y)
y = mean(mean((X-mean(mean(X))).*(Y-mean(mean(Y)))));%/sqrt(mean(mean((X-mean(mean(X))).^2))/sqrt(mean(mean((Y-mean(mean(Y))).^2))));
end

%% Training

function trainPDM(mult,N,alg)
% Training function. Reads the new data format.

% mult = false; % read multiple files
% N=52; % number of points in the model

global trainPDMpathName
cellArray=[];
if ~exist('trainPDMpathName','var'), trainPDMpathName = ''; end
while true
  % Opening a mesh file
  [FileName,PathName] = uigetfile('*.mat','Select File with Mesh Data...','MultiSelect','on',[trainPDMpathName '/']);
  if isequal(FileName,0)
    if size(cellArray,1)<1; 
      return;
    else
      break;
    end
  end
  trainPDMpathName = PathName;
  if ~iscell(FileName), FileName = {FileName}; end
  for u=1:length(FileName)
    load(fullfile2(PathName,FileName{u}),'cellList');

    % Now building an array (cellArray) of the training set points
    for frame=1:length(cellList)
      for cell=1:length(cellList{frame})
        %Convert polygons to contour
        if isempty(cellList{frame}{cell}), continue; end
        mesh=cellList{frame}{cell}.mesh;
        if length(mesh)<4, continue; end
        ctr = [mesh(1:end-1,1:2);flipud(mesh(:,3:4))];
        % ctr = [reshape(plg(1,:,1:end-1),2,[])';plg(2,:,end);plg(1,:,end);plg(3,:,end);flipud(reshape(plg(4,:,1:end-1),2,[])')];
        dctr=diff(ctr,1,1);
        len=cumsum([0;sqrt((dctr.*dctr)*[1;1])]);
        l=length(ctr)-1;
        len1=linspace(0,len(l/2+1),N/2+1);
        len2=linspace(len(l/2+1),len(end),N/2+1);
        len3=[len1(1:end-1) len2];
        ctr1=interp1(len,ctr,len3);
        ctr2=ctr1(2:end,:);%The first and last points are no more the same. The end points are N/2 and N
        ctr2(:,1)=mean(ctr2(:,1))-ctr2(:,1);
        ctr2(:,2)=mean(ctr2(:,2))-ctr2(:,2);
        if len3(end)<5*sqrt(sum((ctr2(N/2,:)-ctr2(N,:)).^2,2))
          cellArray=cat(3,cellArray,ctr2);
        end
      end
    end
  end
  if ~mult, break; end
end
gdisp('Mesh data loaded')
ncells = size(cellArray,3);

% Prealigning the set
time1 = clock;
for i=1:ncells;
  cCell = cellArray(:,:,i);
  alpha = angle(cCell(N,1)+j*cCell(N,2)-(cCell(N/2,1)+j*cCell(N/2,2)));
  cCell = M(cCell,alpha);
  cen = (cCell(ceil(N/4),1)+cCell(ceil(N*3/4),1)+j*cCell(ceil(N/4),2)+j*cCell(ceil(N*3/4),2))/2;
  alpha = angle(cCell(N/2,1)+j*cCell(N/2,2)-cen);
  if alpha>0, cCell(:,2)=-cCell(:,2); cCell=flipud(circshift(cCell,1)); end
  cellArray(:,:,i) = cCell;
end
gdisp('Cells prealigned')

if alg==2

  cellArray2 = cellArray;
  cellArray2(:,2,:) = -cellArray2(:,2,:);
  cellArray2 = flipdim(cellArray2([N 1:N-1],:,:),1);

  cellArray = cat(3,cellArray,cellArray2);
  w=ones(N,1);
  %w2=repmat(w,[1 1 ncells]);
  mCell = cellArray(:,:,1);
  for i=1:10%:10
    %dist=sum(sum((s1-s2).^2,2).*w2,1);
    for k=1:ncells
      cCell = cellArray(:,:,k);
      tmin=fminbnd(@distM,-pi/5,pi/5);
      cellArray(:,:,k) = M(cCell,tmin);
    end
    mCell = mean(cellArray(:,:,:),3);
    w = 1./var(sum(cCell-mCell,2));
    gdisp(['Aligning: Step ' num2str(i)])
  end
  gdisp('Cells aligned')

  % principal components analysis
  data = [reshape(cellArray(:,1,:),N,[]);reshape(cellArray(:,2,:),N,[])]';
  [coefPCA,scorePCA,latPCA] = princomp(data);
  gdisp('PCA completed')

elseif  alg==3
  
  w=ones(N,1);
  %w2=repmat(w,[1 1 ncells]);
    cellArray1 = cellArray;
    cellArray2 = cellArray;
    cellArray2(:,2,:) = -cellArray2(:,2,:);
    cellArray2 = flipdim(reshape(circshift(cellArray2,1),size(cellArray)),1);
    cellArray2 = cat(3,cellArray,cellArray2);
    mCell = mean(cellArray2(:,:,:),3);
    cellParamArray = zeros(length(cellArray),5);
    cellParamArray(:,5)=1;
    % cellArray(:,:,1);
  for i=1:5
    %dist=sum(sum((s1-s2).^2,2).*w2,1);
    for k=1:ncells
      cCell = cellArray(:,:,k);
      cCell(:,1) = cCell(:,1)-mean(cCell(round([N/4 N*3/4]),1));
      cCell(:,2) = cCell(:,2)-mean(cCell(round([N/4 N*3/4]),2));
      tmin = fminsearch(@distParam2cell,cellParamArray(k,3:5));
      cellParamArray(k,3:5) = tmin;
      %cCell = param2cell(cellParamArray(k,3:5))+repmat(cellParamArray(k,1:2),size(cCell,1),1);
      % tmin=fminbnd(@distM,-pi/5,pi/5);
      % cCell = M(cCell,tmin);
      % bmin=fminbnd(@distB,-0.1,0.1);
      % cCell = Bu(cCell,bmin);
      % smin=fminbnd(@distS,0.1,10);
      % cCell = St(cCell,smin);
      cellArray1(:,:,k) = cCell;
    end
    cellArray2 = cellArray1;
    cellArray2(:,2,:) = -cellArray2(:,2,:);
    cellArray2 = flipdim(reshape(circshift(cellArray2,1),size(cellArray)),1);
    cellArray2 = cat(3,cellArray1,cellArray2);
    mCell = mean(cellArray2(:,:,:),3);
    w = 1./var(sum(cCell-mCell,2));
    time2 = clock;
    gdisp(['Aligning: Step ' num2str(i) ', elapsed time ' num2str(etime(time2,time1)) ' s']); 
  end
  gdisp('Cells aligned')

  % principal components analysis
  data = [reshape(cellArray2(:,1,:),N,[]);reshape(cellArray2(:,2,:),N,[])]';
  [coefPCA,scorePCA,latPCA] = princomp(data);
  gdisp('PCA completed')
  
end

% Saving / dislaying data
[FileName,PathName] = uiputfile('*.mat','Select File for PCA Data...');
save(fullfile2(PathName,FileName),'coefPCA','scorePCA','latPCA','mCell');

figure;
plot(mCell(:,1),mCell(:,2),'-',reshape(cellArray2([N/4 N/2 N*3/4 N],1,:),4,[]),reshape(cellArray2([N/4 N/2 N*3/4 N],2,:),4,[]),'.')

function res=distParam2cell(params) % all parameters but shift
  res=sum(sum((mCell-param2cell(params)).^2,2).*w);
  end

function res=param2cell(params) % all parameters but shift
  res=St(Bu(M(cCell,-params(1)),params(2)),params(3));
  end

function res=distM(t) % rotating
  res=sum(sum((M(cCell,t)-mCell).^2,2).*w);
  end

function res=distB(t) % unbending
  res=sum(sum((Bu(cCell,t)-mCell).^2,2).*w);
end

function res=distS(t) % stretching
  res=sum(sum((St(cCell,t)-mCell).^2,2).*w);
  end
end

function b=Bu(a,t)
  % unbends a set of points (up if t>0) assuming initial radius of curvature 1/t
  if abs(t)<1/10000000, b=a; return; end
  R = sign(t)/abs(t);% sign(t)*max(1/abs(t),1.1*max(-a(:,2)));
  tmp = a(:,1) + j*(a(:,2)+R);
  b(:,1) = -(mod(angle(tmp)+sign(t)*pi/2,2*pi)-pi)*R;
  b(:,2) = sign(t)*abs(tmp)-R;
  % r = sqrt((a(:,2)+R).^2+a(:,1).^2);
  % phi = asin(a(:,1)./r);
  % b(:,1) = phi.*R;
  % b(:,2) = r-R;
end

function b=St(a,t)
  % stretches a cell in x direction
  b=a;
  b(:,1) = a(:,1)/t;
end

% The next 2 functions correct a bug in MATLAB uiputfile/uigetfile functions
function [d,e] = uigetfile2(a,b,c)
  f=true;
  while f
    try
      [d,e] = uigetfile(a,b,c);
      f=false;
    catch
      f=true;
    end
  end
end

function [d,e] = uiputfile2(a,b,c)
  f=true;
  while f
    try
      [d,e] = uiputfile(a,b,c);
      f=false;
    catch
      f=true;
    end
  end
end

function res = slashsplit(str)
  % This function returns the highest level folder name from a path
  split = splitstr('/', str);
  split = splitstr('\', split{end});
  res = split{end};
end

function res = slashsplit2(str)
  % This function truncates the path to a folder leaving the two highest
  % level folder names
  split = splitstr('/', str);
  if length(split)>1, split = [split{end-1} '/' split{end}]; else split = split{end}; end
  split = splitstr('\', split);
  if length(split)>1, split = [split{end-1} '\' split{end}]; else split = split{end}; end
  res = split;
end

function parts = splitstr(divider, str)
% This function splits a string into pieces at every occurrence of
% "divider" and returns the result as a cell array of strings. "divider"
% is not included in the output.
   splitlen = length(divider);
   parts = {};
   while 1
    k = strfind(str, divider);
    if isempty(k)
     parts{end+1} = str;
     break
    end
    parts{end+1} = str(1 : k(1)-1);
    str = str(k(1)+splitlen : end);
   end
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

function s = ifdescp(z, nd)
  % S = IFDESCP(Z, ND) computes the inverse Fourier descriptors of Z, which
  % is a sequence of Fourier descriptors arrived at by using FDESCP.
  %
  % "Digital Image Processing Using Matlab", Gonzalez, Woods, Eddins
  % Needed for Fourier smoothing

  np = length(z);
  if nargin==1 || nd>np
    nd = np;
  end
  x = 0:(np-1);
  m = ((-1).^x)';
  d = round((np - nd)/2);
  z(1:d) = 0;
  z(np-d+1:np) = 0;
  zz = ifft(z);
  s(:,1) = real(zz);
  s(:,2) = imag(zz);
  s(:,1) = m.*s(:,1);
  s(:,2) = m.*s(:,2);
end

function z = frdescp(s)
  % Z = FRDESCP(S) computes the Fourier descriptors of S, which is an np-by-2
  % sequence of image coordinates describing a boundary.
  %
  % "Digital Image Processing Using Matlab", Gonzalez, Woods, Eddins.
  % Needed for Fourier smoothing

  [np, nc] = size(s);
  if nc ~=2
    error('S must be of size np-by-2.');
  end
  if np/2 ~= round(np/2);
    s(end+1,:) = s(end, :);
    np = np + 1;
  end

  x = 0:(np-1);
  m = ((-1).^x)';
  s(:,1) = m .* s(:,1);
  s(:,2) = m .* s(:,2);
  s = s(:,1) + sqrt(-1)*s(:,2);
  z = fft(s);
end

function c=circshift(a,b)
  L=size(a,1);
  c = a(mod((0:L-1)-b,L)+1,:);
end

%function gdisp(data)
%  % text display function
%  % alternative to text display to screen
%  % modified from stand-alone version by removing 'CLOSE','VISIBLE','HIDE' commands 
%
%  global gdisphandles logcheckw
%  
%  if logcheckw
%    disp(data)
%  end
%  
%  if ~exist('gdisphandles','var') || ~isfield(gdisphandles,'gdispfig') || isempty(gdisphandles.gdispfig) || ~ishandle(gdisphandles.gdispfig)
%    return
%  end
%  
%  maxlines = 200;
%
%  if ~isa(data,'char'), return; end
%
%  gdisphandles.text = strvcat(gdisphandles.text,data);
%  nlines = size(gdisphandles.text,1);
%  if nlines>maxlines
%    gdisphandles.text = gdisphandles.text(nlines-maxlines+1:nlines,:);
%  end
%  set(gdisphandles.wnd,'String',gdisphandles.text);
%  refresh(gdisphandles.gdispfig)
%  pause(0.005);
%  try
%    gdisphandles.gdispobj.setCaretPosition(gdisphandles.gdispobj.getDocument.getLength);
%  catch
%  end
%end

function gdisp(data)
  % text display function
  % alternative to text display to screen
  % modified from stand-alone version by removing 'CLOSE','VISIBLE','HIDE' commands 
  % SH 120815: Modified to give text output
  
  disp(data);
end



function res = spsmooth(x,y,p,xx)
  % A simple smoothing cubic spline routine. It takes original points y, 
  % their parameterization x, tolerance p, and the new parameterization xx.
  xi = reshape(x,[],1); % make x vertical
  yi = y'; % make y vertival
  n = size(xi,1);
  ny = size(yi,2);
  nn = ones(1,ny);
  dx = diff(xi);
  drv = diff(yi)./dx(:,nn);
  % adx = abs(dx);
  % w = max([adx;0],[0;adx])/mean(adx);
  w = ones(length(x),1);
  if n>2
     idx = 1./dx;
     R = spdiags([dx(2:n-1), 2*(dx(2:n-1)+dx(1:n-2)), dx(1:n-2)], -1:1, n-2, n-2);
     Qt = spdiags([idx(1:n-2), -(idx(2:n-1)+idx(1:n-2)), idx(2:n-1)], 0:2, n-2, n);
     W = spdiags(w,0,n,n);
     Qtw = Qt*spdiags(sqrt(w),0,n,n);
     u = ((6*(1-p))*(Qtw*Qtw.')+p*R)\diff(drv);
     yi = yi - (6*(1-p))*W*diff([zeros(1,ny); diff([zeros(1,ny);u;zeros(1,ny)])./dx(:,nn); zeros(1,ny)]);
     c3 = [zeros(1,ny);p*u;zeros(1,ny)];
     c2 = diff(yi)./dx(:,nn)-dx(:,nn).*(2*c3(1:n-1,:)+c3(2:n,:));
     coefs = reshape([(diff(c3)./dx(:,nn)).',3*c3(1:n-1,:).',c2.',yi(1:n-1,:).'],(n-1)*ny,4);
  else % straight line output
     coefs = [drv.' yi(1,:).'];
  end
  breaks = xi.';
  sizec = size(coefs);
  k = sizec(end);
  l = prod(sizec(1:end-1))/ny;
  [mx,nx] = size(xx);
  lx = mx*nx;
  xs = reshape(xx,1,lx);
  [tmp,index] = histc(xs,[-inf,breaks(2:end-1),inf]);
  NaNx = find(index==0);
  index = min(index,numel(breaks)-1);
  index(NaNx) = 1;
  xs = xs-breaks(index);
  xs = reshape(repmat(xs,ny,1),1,ny*lx);
  index = reshape(repmat(1+ny*index,ny,1)+repmat((-ny:-1).',1,lx), ny*lx, 1 );
  v = coefs(index,1).';
  for i=2:k
     v = xs.*v + coefs(index,i).';
  end
  if ~isempty(NaNx) && k==1 && l>1, v = reshape(v,ny,lx); v(:,NaNx) = NaN; end
  v = reshape(v,ny*mx,nx);
  
  res = reshape(v,[ny,length(xx)]);
end

function res = fullfile2(varargin)
  % This function replaces standard fullfile function in order to correct
  % a MATLAB bug that appears under Mac OS X
  % It produces results identical to fullfile under any other OS
  arg = '';
  for i=1:length(varargin)
    if ~strcmp(varargin{i},'\') && ~strcmp(varargin{i},'/')
      if i>1, arg = [arg ',']; end
      arg = [arg '''' varargin{i} ''''];
    end
  end
  eval(['res = fullfile(' arg ');']);
end

function res = checkparam(p,varargin)
  % this function checks if at least one of the provided parameters is
  % missing
  res = false;
  if isempty(p)
    res = true;
  end
  for i=1:length(varargin)
    if ~isfield(p,varargin{i})
      res = true;
    end
  end
end

function res = readtextfile(filename)
  fid = fopen(filename);
  str = fscanf(fid,'%c');
  fclose(fid);
  newline = [0 regexp(str,'\n') length(str)];
  res = {};
  for i=2:length(newline)
    res = [res str(newline(i-1)+1:newline(i)-1)];
  end
end

function initP
  global se maskdx maskdy
  
  % Global derived parameters
  % initModel
  
  se = strel('arb',[0 1 0;1 1 1;0 1 0]); % erosion mask, can be 4 or 8 neighbors
  maskdx = fliplr(fspecial('sobel')'); %[-1 0 1; -2 0 2; -1 0 1]/2; % masks for computing x & y derivatives
  maskdy = fspecial('sobel');%[1 2 1; 0 0 0; -1 -2 -1]/2;
end


function initMTrack

  % define/redefine globals
  global rmask cellList cellListN selectedList shiftframes imsizes imageFolders imageLimits handles logcheckw regionSelectionRect shiftframes shiftfluo%#ok<REDEF>
  initP;
  cellList = {{}}; % list to store fitted cells
  cellListN = [];
  selectedList = [];
  imsizes = zeros(5,3);
  regionSelectionRect = [];
  shiftframes = [];
  imageFolders = {'','','',''};
  logcheckw = true;
  imageLimits = {[0 1],[0 1],[0 1],[0 1]};
  shiftfluo = [0 0; 0 0];
  
  rmask{1} = strel('arbitrary',[0 0 0 0 0; 0 0 0 0 0; 1 1 0 1 1; 0 0 0 0 0; 0 0 0 0 0]);
  rmask{2} = strel('arbitrary',[1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 1 0; 0 0 0 0 1]);
  rmask{3} = strel('arbitrary',[0 0 1 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 1 0 0; 0 0 1 0 0]);
  rmask{4} = strel('arbitrary',[0 0 0 0 1; 0 0 0 1 0; 0 0 0 0 0; 0 1 0 0 0; 1 0 0 0 0]);
   
end


%---------------------------------------------------
function position = assignSegment(x,y,mesh,box,coordOffset);
  %function position = assignSegment(x,y,mesh,box);
  % Assign segment number in which spot is located inside cell
  % Written: SH120817
  if ~exist('box','var')
    box = [1 1];
  end
  if ~exist('coordOffset','var')
    coordOffset= 1;
  end

  
  position = zeros(size(x));
  
  for kk=1:size(mesh,1)-1
    plgx = [mesh(kk,[1 3]) mesh(kk+1,[3 1])]-box(1)+1;
    plgy = [mesh(kk,[2 4]) mesh(kk+1,[4 2])]-box(2)+1;
    inPol = inpolygon(x,y,plgx,plgy);
    inPol = inPol*kk;
    position = position + inPol;
  
  end
end


%----------------------------------------------
function [XIn] = dblTformInv(tformDbl,XBase)

  xBase = XBase(:,1);
  yBase = XBase(:,2);
  [x1,y1] = tforminv(tformDbl{1},xBase,yBase);
  [xIn,yIn] = tforminv(tformDbl{2},x1,y1);
  XIn =[xIn,yIn];

end
%-----------------------------------------------
function dblTform_InToBase = getDblTform(In,Base,TFORM_TYPE) 

  dblTform_InToBase{1} = cp2tform(In,Base,'similarity');
  [XBase1,YBase1] =tforminv(dblTform_InToBase{1},Base(:,1),Base(:,2));
  %dblTform_InToBase{2} = cp2tform(In,[XBase1,YBase1],'similarity');
  %dblTform_InToBase{2} = cp2tform(In,[XBase1,YBase1],'piecewise linear');
  dblTform_InToBase{2} = cp2tform(In,[XBase1,YBase1],TFORM_TYPE);

end

