function [chk,h1,h2] = hb_nii_verify_space_match(f1,f2)
% HB_NII_VERIFY_SPACE_MATCH verifies match between two nifti files in terms
% of their 3D dimension and coordinate space. That is, whether the two
% vloumes are matched voxel-to-voxel in space. The two inputs can be either
% file paths or nifti headers.
%
% Note: If f1, f2, or both are 4D, the assumption is that all frames have a
% matching header, and thus, only the first frame in the AD volume will be
% checked.
% 
% Hamid Behjat

h1 = gethead(f1);
h2 = gethead(f2);
chk1 = isequal(h1.dim, h2.dim);
chk2 = all(abs(h1.mat-h2.mat)<1e-6,'all');
chk = chk1 && chk2;

end



%==========================================================================
function h = gethead(f)
if ischar(f)
    if exist(f,'file')
        if endsWith(f, '.nii.gz')
            gunzip(f);
            d = strrep(f, '.gz', '');
            h = spm_vol(d);
            delete(d);
        elseif endsWith(f, '.nii')
            h = spm_vol(f);
        else
            error('Incorrect file format.');
        end
    else
        fgz = [f,'.gz'];
        if exist(fgz, 'file')
            gunzip(fgz);
            h = spm_vol(f);
            delete(f);
        else
            error('Missing file: %s',f);
        end
    end
else
    assert(isstruct(f));
    h = f;
end
if length(h)>1
    h = h(1);
end
end