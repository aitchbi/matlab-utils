function [r,s,I,sts] = hb_get_voxres(m)
% HB_GET_VOXRES determines voxels resolution from nifti mat file. 
% m: h.mat of nifti header; 4x4 matrix
% r: voxel size
% s: sign aociated to each element of r in m
% I: 3x3 matrix showing position of non-zero elements m(1:3,1:3)
% sts: logical; is voxel isotropic
%
% Hamid Behjat

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
    sts = 'Voxels are non-isotropic.';
else
    sts = 1; 
end
end