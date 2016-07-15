function limVal= percentileLim(data,tol)
% find the nth% of the data
% tol takes values between 0,1
% limVal is the nearest point the the dataset to the specified limits

%remove nans
data(isnan(data)) =[];

data = sort(data(:));
limIdx = max(1,round(tol*numel(data)));
limVal= data(limIdx);

