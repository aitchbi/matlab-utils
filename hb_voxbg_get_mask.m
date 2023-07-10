function f_mask = hb_voxbg_get_mask(f_graph)
% HB_VOXBG_GET_MASK extracts 3D graph mask from a voxBG and writes it as a
% nifti file in teh same directory as the graph.
%
% Inputs:
%   f_g: voxel-wise brain graph; .mat file. 
%
%  Outputs:
%   f_m: nifit binary mask file, showing indices in 3D volume corresponding
%   to graph vertices. The file is written in the same folder as f_g, and
%   the name will be the same as f_g appended with '_mask'.
%
% Examples:
% f = hb_voxbg_get_mask('/a/b/c/G.gmlh.res1250.spaceT1w.mat');
%
% Dependencies: 
% .SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12
%
% Hamid Behjat

assert(endsWith(f_graph,'.mat'));

f_mask = strrep(f_graph,'.mat','_mask.nii');

d = load(f_graph);
G = d.G;

h = struct;
h.dim = G.dim;
h.mat = G.mat;
h.dt = [16 1];
h.fname = f_mask;
v = zeros(h.dim);
v(G.indices) = 1;
spm_write_vol(h,v);

gzip(f_mask);
delete(f_mask);

f_mask = strcat(f_mask,'.gz');

end