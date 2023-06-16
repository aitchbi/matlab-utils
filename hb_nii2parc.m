function [P, l] = hb_nii2parc(f_i, f_p, varargin)
% HB_NII2PARC parcellates a given nifti based on a given parcellation.
%
% 
% Inputs:
%   f_i: a nifti volume; current code only works for 3D volume, not 4D.
% 
%   f_p: a segmentation/parcellation nifti file, e.g. a mask or a
%   cortical parcellation.
%
% Outputs:
%   P: an MxF cell array, wherein array (i,j) stores voxel values of j-th
%   frame of f_i associated to label l(j) in f_p. The array is sorted based
%   on sorted labels in f_p; i.e., P{1} and P{end} contain voxels in f_i
%   associated to the smallest and largest non-zero label in f_p,
%   respectively.
%
%   l: the labels in f_p associated to each element of P; the labels are
%   sorted in ascending order, i.e, L(i) < L(i+1) for all i \in [1, Np]
%   where Np denotes the number of unique non-zero labels in f_p.
%
% 
% Dependencies:
% .SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_nii_reslice.m 
% .https://github.com/aitchbi/matlab-utils/blob/main/spm_modified
% 
% HB

sts = chkmatch(f_i, f_p);
assert(length(spm_vol(f_i))==1, 'extend for 4D');

% load f_i
if sts
    % f_i & f_p match voxel-to-voxel
    v_i = spm_read_vols(spm_vol(f_i));
else
    % reslice f_i to match f_p voxel-to-voxel
    f_tmp = strrep(f_i, '.nii', '_tmp.nii');
    hb_nii_reslice(f_i,f_p,[],f_tmp);
    assert(chkmatch(f_tmp, f_p));
    v_i = spm_read_vols(spm_vol(f_tmp));
    delete(f_tmp);
end

v_p = spm_read_vols(spm_vol(f_p));
l = sort(unique(v_p));
assert(l(1)==0);
l = l(2:end); % exclude 0
Np = length(l); 
P = cell(Np,1);
for iP=1:Np
    P{iP} = v_i(v_p==l(iP));
end
end

%==========================================================================
function sts = chkmatch(f1,f2)
tol = 1e-6;
h1 = spm_vol(f1);
h2 = spm_vol(f2);
if length(h1)>1
    h1 = h1(1);
end
if length(h2)>1
    h2 = h2(1);
end
d1 = isequal(h1.dim, h2.dim);
d2 = all(all(abs(h1.mat-h2.mat)<tol));
if d1 && d2
    sts = true;
else
    sts = false;
end
end