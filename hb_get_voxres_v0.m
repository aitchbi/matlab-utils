% [07.12.2025] this version cannot handle tilted unit vectors; ie it throws
% an error if the direction vectors are not aligned with x, y, and z axes
% (ie each column in I has more than one non-zero element). for the a
% genralized mat matrix, see hb_get_voxres.m

function [r,s,I,sts] = hb_get_voxres_v0(m)
% HB_GET_VOXRES determines voxels resolution from nifti mat file. 
% m: h.mat of nifti header; 4x4 matrix
% r: voxel size
% s: sign asociated to each element of r in m
% I: 3x3 matrix showing position of non-zero elements m(1:3,1:3)
% sts: logical; is voxel isotropic?
%
% h behjat

m3 = m(1:3,1:3);
I = abs(m3-0)>1e-3;
assert(nnz(I)==3,'mat file has off diagonal elements.');
r = zeros(3,1);
r(1) = m3(1,I(1,:)); % x
r(2) = m3(2,I(2,:)); % y
r(3) = m3(3,I(3,:)); % z
s = sign(r);
r = abs(r); 
d = abs(r-r(1))<1e-3;
if not(all(d))
    sts = 'voxels are non-isotropic.';
else
    sts = 1; 
end
end