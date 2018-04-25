function output=importRapidStorm(fname)
% output=importRapidStorm(fname) is a function to import a rapidstorm
% localization file.
% (Tried using the import function from PALMsiever, but it had java errors)

% input check
if ~ischar(fname)
    error('input should be a string')
end

% open file
fid=fopen(fname);
if fid == 0
    error('failed to loadfile')
end

% prealloc output?
output=cell(1,2);

% read lines
lineNr=0;
while feof(fid)== 0
    lineNr=lineNr+1;
    aline=fgetl(fid);
    output(1,lineNr)=aline;
end

% close the file
closed=fclose(fid);
if closed ~= 0
    disp('failed to close file')
end
    
end