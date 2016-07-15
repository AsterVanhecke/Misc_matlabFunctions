function summaryDataStruct= readImSummaryMM(fname)

fid = fopen(fname,'r');
%skip pre summary part of header
skipN = 40;
for i=1:skipN
  A = fscanf(fid,'%c',1);
end

%preallocate long string
nPre=1e4;
summaryData = cast(zeros(1,nPre),'char');

endOfSumData = false;
ii =1;
while endOfSumData ==false
  summaryData(ii) = fscanf(fid,'%c',1); 
  if summaryData(ii) == '}'
    endOfSumData = true;
  end
  ii = ii + 1;
end
summaryData(ii:end)=[];

summaryDataStruct = loadjson(summaryData);
