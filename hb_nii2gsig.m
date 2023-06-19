function S = hb_nii2gsig(f, I)
%
% Inputs:
%   f: full path of nifti file
%   I: indices in 3D volume 
% Outputs: 
%   S: graph signals, one per column.
%
% Dependencies: 
% . SPM12 package
% . hb_nii_load.m
%
% HB

if endsWith(f, '.gz')
    fgz = f;
    f = strrep(f, '.nii.gz', '.nii');
else
    fgz = [f, '.gz'];
end

if exist(f, 'file')
    CleanUp = false;
else
    if exist(fgz, 'file')
        gunzip(fgz);
        CleanUp = true;
    else
        error('file not found');
    end
end

[v, h] = hb_nii_load(f, 'IndicesToLoad', I);

Ng = length(I);

Nv = length(h);

S = zeros(Ng, Nv);
for k=1:Nv
    d = v(:,:,:,k);
    S(:, k) = d(I) ;
end

if CleanUp
    delete(f);
end
end