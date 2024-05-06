function [volume, info] = hb_get_surfvol(f_surf, f_brain, f_out, varargin)
% HB_GET_SURFVOL extract volume encapsulated by a closed cortical surface. 
%
% Inputs: 
%   f_surf: pathname to surface file (gifti); lh/rh.pial/white.surf.gii
%   f_brain: pathname to subject's skull-stripped brain (nifti)
%   f_out: pathname to output file to write
%
% Outputs: 
%   volume: volume of extracted surface in mm cubic.
%   info: info about the procedure e.g. to visualise detached mask.  
%
% Dependencies: 
% -github.com/aitchbi/matlab-utils [*]
% -github.com/aitchbi/GRASS/tree/main/utils/prune [**]
% -https://www.fil.ion.ucl.ac.uk/spm/software/spm12 [*]
%
% [*] download then addpath(.)
% [**] download then addpath(genpath(.))
%
% Hamid Behjat 

d = inputParser;
addParameter(d,'Parallelize', true);
addParameter(d,'OverwriteExisting', false);
parse(d,varargin{:});
opts = d.Results;

if ~opts.OverwriteExisting
    if exist(f_out,'file') || exist([f_out,'.gz'],'file')
        return;
    end
end

h_mask = spm_vol(f_brain);
mask = spm_read_vols(h_mask);
mask = logical(mask);
mask = hb_make_connected(mask,26);
indices = find(mask);

% adjacencies across whole brain mask
A0 = hb_get_adjacency(mask,26);

% remove edges passing through white surface
if ~iscell(f_surf)
    f_surf = {f_surf};
end
surface = ml_batch_gifti(f_surf,h_mask.mat);
A = ml_prune_adjacency(A0,mask,surface,'parallelize',opts.Parallelize);

% get connected comps
bins = conncomp(graph(A));

% sort comps based on size
bin_count = max(bins);
bins_size = histcounts(bins,1:bin_count+1);

[~,d] = sort(bins_size,'descend');
d = ismember(bins,d(2));
inds_rm = indices(~d); % voxels to remove
mask_marked = mask;  
mask_marked(inds_rm) = 2; % volume showing removed voxels
mask(inds_rm) = 0;

%-Make mask connected.
d1 = nnz(mask);
mask = hb_make_connected(mask,26);
d2 = nnz(mask);
assert(d2>.95*d1,'fishy; too many voxels removed; check resulting mask');

%-Write mask.
h = struct();
h.fname = f_out;
h.dim = h_mask.dim;
h.dt = [spm_type('uint8') 0];
h.mat = h_mask.mat;
spm_write_vol(h,mask);
voxnum = nnz(mask);
voxvol = getvoxvol(h_mask.mat);
volume = voxnum*voxvol;

%-Write base mask with removed parts marked.
[p,n,e] = fileparts(f_out);
h.fname = fullfile(p,[n,'.base_mask_with_removed_part_marked',e]);
spm_write_vol(h,mask_marked);

gzip(h.fname);
delete(h.fname);

%-Write info file.
info = struct;
info.pruning.A_diff_pre_post_white_pruning = A0-A;
info.pruning.mask_with_removed_part_marked = h.fname;
info.pruning.info = fullfile(p,[n,'.extracting_white_structure_info.txt']);
fid = fopen(info.pruning.info,'wt');
fprintf(fid,'=========================================================\n');
fprintf(fid,'%%To plot prunning done to detach surface mask: \n');
fprintf(fid,'%%A0: pre-pruning A. \n');
fprintf(fid,'%%Ap : post-pruning A. \n');
fprintf(fid,'\n');
fprintf(fid,'d = hb_gunzip(info.pruning.mask_with_removed_part_marked); \n');
fprintf(fid,'mask = spm_read_vols(spm_vol(d)); \n');
fprintf(fid,'A0 = hb_get_adjacency(logical(mask),26); \n');
fprintf(fid,'Adiff = info.pruning.A_diff_pre_post_white_pruning; \n');
fprintf(fid,'Ap = A0-Adiff;\n');
fprintf(fid,'sl = 88;       %%axial slice number \n');
fprintf(fid,'d1 = ''grid''; %%''diff'',''none'' \n');
fprintf(fid,'d2 = true;     %%false \n');
fprintf(fid,'d3 = ''r'';    %%color \n');
fprintf(fid,'figure; \n');
fprintf(fid,'hm_plot_adjacency_diff(A0,Ap,mask,sl,d1,d2,d3);\n');
fprintf(fid,'%%or \n');
fprintf(fid,'figure; \n');
fprintf(fid,'hm_plot_adjacency_diff(A0,Ap,logical(mask),sl,d1,d2,d3);\n');
fprintf(fid,'=========================================================\n');
fclose(fid);
end

%==========================================================================
function v = getvoxvol(hmat)
r = hb_get_voxres(hmat);
v = prod(r); % voxel volume in mm^3
end
