function [v_o,n,f,v_rm] = hb_make_connected(v_i,conn)
%HB_MAKE_CONNECTED makes input mask (bw image. 2D or 3D) connected.
% The largest component will be kept and the remainder of voxels/components  
% will be deleted. 
% 
% Inputs:
%   v_i: bw 2D image or 3D volume. Or spm_vol header.
%   conn: 4 or 8 for 2D. 6, 18 or 26 for 3D. See bwconncomp.m for details.
%
% Outputs:
%   v_o: bw image/volume made single connected. 
%   n: number of voxels removed to make mask single connected.
%   f: fraction of voxels removed to make mask single connected.
%   v_rm: bw image/volume showing removed pixels/voxels. 
%
% Dependencies: SPM12 toolbox.
%
% Hamid Behjat 

if isstruct(v_i)
    v_i = spm_read_vols(v_i);
end

vdim = size(v_i);

CC = bwconncomp(v_i,conn);

[~,imax] = max(cellfun(@length,CC.PixelIdxList));

v_o = zeros(vdim);
v_o(CC.PixelIdxList{imax}) = 1;

v_rm = zeros(vdim);
for n=1:length(CC.PixelIdxList)
    if n==imax
        continue;
    end
    v_rm(CC.PixelIdxList{n}) = 1;
end

n = nnz(v_rm); %nnz(v_i)-nmax; 
f = n/nnz(v_i);  
if f>0.05
    warning('A large number of voxels were removed!'); 
end
end