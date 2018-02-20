function VECLIST = readSPNVecList(filename)

fid = fopen(filename);

inputLine = fgetl(fid);
if ~strcmp(inputLine, 'SPN_VECLIST:txt')
    error('first line has to be ''SPN_VECLIST:txt''')
end

inputLine = fgetl(fid);
numVecs = str2double(inputLine);
if isnan(numVecs)
    error('malformatted file')
end

VECLIST = cell(numVecs,1);

for k=1:numVecs
    inputLine = fgetl(fid);
    curVec = str2num(inputLine);
    VECLIST{k} = curVec;
end

fclose(fid);
