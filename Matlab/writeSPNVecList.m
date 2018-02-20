function writeSPNVecList(filename, VECLIST)

if iscell(VECLIST)
    for k=1:size(VECLIST,1)
        if ~isvector(VECLIST{k})
            error('VECLIST must contain only vectors')
        end
    end
end


fid = fopen(filename, 'w');

fprintf(fid, 'SPN_VECLIST:txt\n');

if iscell(VECLIST)
    fprintf(fid, '%d\n', length(VECLIST));
    
    for k=1:numel(VECLIST)
        fprintf(fid, '%d ', VECLIST{k});
        fprintf(fid, '\n');
    end
else
    fprintf(fid, '%d\n', size(VECLIST,1));
    
    for k=1:size(VECLIST,1)
        fprintf(fid, '%d ', VECLIST(k,:));
        fprintf(fid, '\n');
    end
end

fclose(fid);
