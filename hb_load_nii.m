function [v,h] = hb_load_nii(f,headertype)
%HB_LOAD_NII load nifti volume and header using external software. 

if ~exist('headertype','var')
    headertype = 'spm';
end

if contains(f,'.gz')
    gunzip(f);
    f = f(1:end-3);
    DELNONZIP = true;
else
    if exist(f,'file')
        DELNONZIP = false;
    else
        fgz = [f,'.gz'];
        if exist(fgz,'file')
            gunzip(fgz);
            DELNONZIP = true;
        else
            error('File missing: %s',f);
        end
    end
end

switch headertype
    case 'spm'
        assert(exist('spm_vol.m','file'),'"spm" package not in path.')
        h = spm_vol(f);
        v = spm_read_vols(h);
        if DELNONZIP
            delete(f);
        end
    otherwise
        error('extend.')
end
end