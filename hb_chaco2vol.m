function [inds,CC] = hb_chaco2vol(f_src, f_clustTxt, N, f_clustVol, indices, chacoSettings, underlay, clustVolume)
% HB_CHACO2VOL transforms .txt output from Chaco to a volume and saves it
% as a nifti file.
%
% Inputs:
%   f_src: full address to source file (mask) that  has been clustered.
%   [source file should be a single connected binary mask]
%
%   f_clustTxt: This is the output file from Chaco. The value on each row 
%   i specifies which cluster graph vertex i belongs, which is a value in 
%   the range 0:N-1. 
%
%   N: Number of clsuters that Chaco has clustered the source into. 
% 
%   f_clustVol: Full address of the output nifti volume which will be save 
%   on disc. Each clustered region is marked with a value. 
%
%   indices: graph indices in f_src. 
%
%   chacoSettings; (opt) for description.
% 
%   underlay: (opt) ovelay clusters on an underlay and save nii.
%
%   clustVolume: (opt) volume of each cluster in mm^3; for description.
%
% 
% Dependencies: 
%   SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12)
%
% Hamid Behjat

if ~exist('chacoSettings','var')
    chacoSettings = [];
end
if ~exist('underlay','var')
    underlay = [];
end
if ~exist('clustVolume','var')
    clustVolume = [];
end

%-Load cluster indices.
T = readtable(f_clustTxt);
C = T.Var1;
C = C+1; % cluster numbers from Chaco are 0:N-1

%-Sanity check.
assert(numel(unique(C))==N,...
    sprintf('Number of clusters in file (%d) not equal to N (%d)',...
    length(unique(C)),N));

%-Modify colormap.
CC = zeros(size(C));
inds = cell(N,1);
rng('default');
rng(1);
colMap = randperm(N);
for i = 1:N
    d = find(C==i);
    inds{i} = d;
    CC(d) = colMap(i);
end

%-Save volume.
h = spm_vol(f_src);
v = spm_read_vols(h);
assert(isequal(indices(:),find(v(:))));
v_clusters = zeros(size(v));
v_clusters(indices) = CC;
if ~isempty(underlay)
    d = spm_read_vols(spm_vol(underlay));
    d = find(d(:));
    d = setdiff(d,indices);
    v_clusters(d) = N+1;
end
h.fname = f_clustVol;
if ~isempty(clustVolume) && ~isempty(chacoSettings)
    h.descrip = ['clustered (Chaco inputs:',...
        num2str(chacoSettings(1)),'-',...
        num2str(chacoSettings(2)),'-',...
        num2str(chacoSettings(3)),'-',...
        num2str(chacoSettings(4)),'-',...
        num2str(chacoSettings(5)),...
        '); vol of each cluster: ',num2str(clustVolume),' mm^3)'];
elseif ~isempty(clustVolume)
    h.descrip = ['clustered using Chaco (',num2str(N),...
        ' clusters) - volume of each cluster: ',...
        num2str(clustVolume),' mm^3)'];
else
    h.descrip = ['clustered using Chaco (',num2str(N),' clusters)'];
end
spm_write_vol(h,v_clusters);
end
