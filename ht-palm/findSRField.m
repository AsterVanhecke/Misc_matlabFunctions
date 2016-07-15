function datacol = findSRField(parInfo,attributeName, value)

datacol = -1;
nField = numel(parInfo.field);
for ii = 1:nField
  atValue = getfield(parInfo.field(ii).ATTRIBUTE,attributeName);
  if strcmp(atValue,value)
    datacol = ii;
  end
end
