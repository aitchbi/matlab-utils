function [r, voxtype] = hb_get_voxres(m)
% HB_GET_VOXRES determines voxels resolution from input nifti mat file. 
% m: h.mat of nifti header; 4x4 matrix
% r: voxel size along x, y, z axes; 1x3 array
%
% h behjat

%-voxel resolution in each direction---------------------------------------
m3 = m(1:3,1:3);

r = vecnorm(m3, 2, 1); 
% 2: p, p-norm
% 1: dim, ie columns 
%
% each column represents a direction vector (in mm) for stepping one voxel
% along an axis:
%
% 1st colum: x-axis
% 2nd colum: y-axis
% 3rd colum: z-axis
%
% the Euclidean norm of each vector gives the voxel size along each axis

%-isotropic voxels?--------------------------------------------------------
chk1 = abs(r(2)-r(1))<1e-2;
chk2 = abs(r(3)-r(1))<1e-2;
chk3 = abs(r(3)-r(2))<1e-2;
if all([chk1, chk2, chk3])
    voxtype = 'isotropic';
else
    voxtype = 'non-isotropic';
end

end