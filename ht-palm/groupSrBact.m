function srInMeshGrouped = groupSrBact(srInMesh, maxBlink, maxDistance);

nCell = numel(srInMesh);
for ii = 1:nCell
  srInMeshGrouped{ii} = groupBact(srInMesh{ii},maxBlink,maxDistance);
end

%-------------------------------------------
function bactGrouped = groupBact(bact,maxBlink,maxDistance);
parInfo = bact.localizations.parInfo;

xCol = findSRField(parInfo,'semantic','position in sample space in x dimension');
yCol = findSRField(parInfo,'semantic','position in sample space in y dimension');
frameCol =  findSRField(parInfo,'semantic','frame number');

locData = bact.localizations.data;
x = locData(:,xCol);
y = locData(:,yCol);
frame = locData(:,frameCol);

groupId = groupData(x,y,frame,maxBlink,maxDistance);
%just take the average of the parameters (eg dont sum photon counts for now)
bactGrouped = bact;
bactGrouped.l = combineGroup(bact.l,groupId);
bactGrouped.d = combineGroup(bact.d,groupId);
bactGrouped.position = combineGroup(bact.position,groupId);
bactGrouped.localizations.data = combineGroup(bact.localizations.data,groupId);
%-------------------------------------------
function group_id = groupData(x,y,frame,maxF,maxD)
%This is from thomas's PALMsiever grouping function
group_id = zeros(size(x,1),1);

maxD2 = maxD^2;

if maxF>max(frame)
    maxF=1;
end
cons = zeros(size(x));
for ff = min(frame):max(frame)-maxF
    ss = frame>=ff & frame<=(ff+maxF);
    
    iss=find(ss)';
    
    xx = x(iss);
    yy = y(iss);
    
    for ii=1:length(iss)-1
        [ mind2 imind2 ]= min((xx(ii)-xx(ii+1:end)).^2+(yy(ii)-yy(ii+1:end)).^2);
        if mind2 < maxD2
            cons(iss(ii))=iss(ii+imind2);
        end
    end
    
end

% Remove self
cons(cons'==(1:length(cons)))=0;

% % Assign ID
% cur = 1; curID = 1; 
% IDs = zeros(N,1);
% while cur <= N
%     if IDs(cur) > 0 % already assigned to a group
%         if cons(cur)==0 % last member of the group
%             IDs(cons(cur))=IDs(cur); 
%         end
%     else
%         IDs(cur)=curID;
%         curID=curID+1;
%     end    
%     cur = cur + 1;
% end

% Calculate group frame IDs
N = size(cons,1);
curfg=1; group_frame_id=ones(N,1); group_id_=1:N;
while sum(group_frame_id==curfg)>0
    group_frame_id(cons(group_frame_id==curfg & cons>0))=group_frame_id(cons(group_frame_id==curfg & cons>0))+1;
    group_id_(cons(group_frame_id==curfg & cons>0))=group_id_(group_frame_id==curfg & cons>0);
    curfg=curfg+1;
end

[zzz_ zzz_ group_id]= unique(group_id_);


%%-------------------------------------------
function dataGrouped = combineGroup(data,groupId);
nCol = size(data,2);
nGroup = numel(unique(groupId));
dataGrouped = zeros(nGroup,nCol);

for ii = 1:nCol
  dataGrouped(:,ii) = groupBy(data(:,ii),groupId(:));
end

%%-------------------------------------------
function gX = groupBy(X, ID)
%Thomas's palmsiever groupBy function
[a a ID] = unique(ID);
gX = accumarray(ID,X)./accumarray(ID,1);
