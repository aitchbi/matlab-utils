function [f_gparc, P, f_nolabel, I, V] = hb_voxbg_get_maskedparc(f_graph,f_parc,varargin)
% HB_VOXBG_GET_MASKEDPARC generates a masked version of input parcellation
% such that only voxels that fall within the graph mask are returned.
%
% Inputs:
%   f_graph: voxel-wise brain graph; .mat file. 
%
%   f_parc : a parcellation of the subjects brain, defined in the same in
%   which the graph was defined. It is assumed that the parcellation file
%   is in register with the space in which the graph was defined. However,
%   if it is does not voxel-to-voxel match with the space in which tha
%   graph was defined (i.e. its dim and mat as given in its header not
%   matching that of the graph), it will be resliced to match the space in
%   which the graph was defined.
%
%   Name-Value Pair Arguments (optional):
%   GzipOutput: logical, specifying whether or not to write output file as
%   gzip; default: true.
%
%   WhichLabels: a vecor listing the lables from f_parc to consider; labels
%   in f_parc that are not included in this list will be excluded, thus, no
%   returned in output file. Default: all labels within f_parc will be
%   consiudered.
%
%   OutputFileName: nifti file name/path for output file. If not given, a
%   name is generated based on f_graph and f_parc, and the file is written
%   in same directory as f_graph. File format: *.nii or *.nii.gz.
%
%   WriteFileShowingGraphMaskVoxelsWithNoLabel: logical, specifying whether
%   or not to write a volume in which voxels within the garph mask that had
%   no label in f_parc are marked; default: false.
% 
%  Outputs:
%   f_gparc: f_parc masked by the graph's mask; i.e., only voxels that are
%   associated with graph vertices are returned as labels. Voxels in the
%   graph mask that did not fall within f_parc are retuned as 0, just like
%   all other voxels that fall outside the graph mask.
%
%   P: list of parcel labels included in generated file, f_gparc. This list
%   is not necessarily equal to the list of labels in the input file
%   f_parc, because: (i) some labels in f_parc might belong to voxels that
%   are not within the graph mask, (ii) if optional input 'WhichLabels'
%   given, potentially only a subset of labels from f_parc are considered. 
%
%   f_nolabel: file showing voxels from graph mask that were not assigned a
%   label. Generetaed if WriteFileShowingGraphMaskVoxelsWithNoLabel = true.
% 
%   I: graph indices* that were not assigned a label since they were not
%   labelled in f_parc. *a vector of intergers from the the set {1:N},
%   where N denotes the size of the graph. If [], all graph vertices (i.e.,
%   voxels in graph's mask) were labelled.
%
%   V: voxel indices in the graph mask that were not assigned a label. This
%   is a vector with teh same length as I. V(k) is the voxel index
%   associated to graph vertex I(k), i.e. V(k) = G.indices(I(k)).
%
% Examples:
% [f_gp, P] = hb_voxbg_get_maskedparc(f_g,f_p);
% [f_gp, P, ~, I, V] = hb_voxbg_get_maskedparc(f_g,f_p);
% [__, f_nolabel] = hb_voxbg_get_maskedparc(__, 'WriteFileShowingGraphMaskVoxelsWithNoLabel', true);
% hb_voxbg_get_maskedparc(__, 'WhichLabels', 1:200);
% hb_voxbg_get_maskedparc(__, 'OutputFileName', 'f.nii');
%
% Dependencies: 
% .SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12
% .https://github.com/aitchbi/matlab-utils/hb_nii_reslice.m
%
% Hamid Behjat

assert(endsWith(f_graph,'.mat'));

d = inputParser;
addParameter(d,'GzipOutput', true);
addParameter(d,'WhichLabels', []);
addParameter(d,'OutputFileName', []);
addParameter(d,'WriteFileShowingGraphMaskVoxelsWithNoLabel', false);
parse(d,varargin{:});
opts = d.Results;

if not(isempty(opts.OutputFileName))
    if endsWith(opts.OutputFileName, '.gz')
        assert(isequal(opts.GzipOutput,true),...
            printf('Discrepency between optional inputs %s and %s.',...
            'OutputFileName',...
            'GzipOutput'));
    end
end

%-Load graph.
d = load(f_graph);
G = d.G;
h_g = struct;
h_g.dim = G.dim;
h_g.mat = G.mat;

%-Build graph mask.
f_mask = strrep(f_graph,'.mat','___tmp_mask.nii');
F = appendcleanup([],f_mask);
h_m = struct;
h_m.dim = h_g.dim;
h_m.mat = h_g.mat;
h_m.dt = [64 0];
h_m.fname = f_mask;
v_m = zeros(h_m.dim);
v_m(G.indices) = 1;
spm_write_vol(h_m,v_m);

%-Load & verify parcellation.
h_p = spm_vol(f_parc);
chk = verifymatch(h_p,h_g);
if chk
    v_p = spm_read_vols(h_p);
else
    % reslice parcellation to match graph mask, voxel-by-voxel.
    if endsWith(f_parc, '.nii')
        f_p = f_parc;
    elseif endsWith(f_parc, '.nii.gz')
        gunzip(f_parc);
        f_p = strrep(f_parc,'.gz','');
        F = appendcleanup(F,f_p);
    end
    f_ptmp = strrep(f_p,'.nii','___tmp.nii');
    hb_nii_reslice(f_p,f_mask,0,f_ptmp);
    h_atmp = spm_vol(f_ptmp);
    assert(verifymatch(h_atmp,h_g));
    F = appendcleanup(F,f_ptmp);
    v_p = spm_read_vols(h_atmp);
end
P = sort(unique(v_p));
assert(P(1)==0);
P = P(2:end);
if not(isempty(opts.WhichLabels))
    d = ismember(opts.WhichLabels,P);
    if not(all(d))
        nd = not(d);
        fprintf('\nMissing labels in f_parc: %d', opts.WhichLabels(nd));
        error('%d specified labels missing in f_parc.',nnz(nd));
    end
    P = opts.WhichLabels;
end
Np = length(P);

%-Parcellation in graph mask.
if isempty(opts.OutputFileName)
    [~,d] = fileparts(f_parc);
    f_gparc = strrep(f_graph, '.mat', sprintf('_%s.nii', d));
else
    f_gparc = opts.OutputFileName;
    if endsWith(f_gparc, '.nii.gz')
        f_gparc = strrep(f_gparc, '.gz', '');
        opts.GzipOutput = true;
    else
        assert(endsWith(f_gparc, '.nii'));
    end
end
v_mp = zeros(size(v_m));
for k=1:Np
    v_mp(and(v_p==P(k),v_m==1)) = P(k);
end
h_mp = struct;
h_mp.fname = f_gparc; 
h_mp.dim   = h_g.dim; 
h_mp.mat   = h_g.mat;
h_mp.dt    = [64 0]; 
spm_write_vol(h_mp,v_mp);
if opts.GzipOutput
    gzip(f_gparc);
    F = appendcleanup(F,f_gparc);
end
Im  = find(v_m);
Imp = find(v_mp);
if isequal(Im,Imp)
    I = [];
    V = [];
else
    V = Im(not(ismember(Im,Imp)));
    [~,I] = ismember(V,G.indices);
    if opts.WriteFileShowingGraphMaskVoxelsWithNoLabel
        [d1,d2] = fileparts(f_gparc);
        if endsWith(d2,'.nii')
            d2 = strrep(d2,'.nii','');
        end
        f_nolabel = fullfile(d1,...
            [d2,'_unlabelled_graph_mask_voxels.nii']);
        h_nl = struct;
        h_nl.fname = f_nolabel;
        h_nl.dim = h_g.dim;
        h_nl.mat = h_g.mat;
        h_nl.dt = [64 0];
        v_nl = zeros(size(v_m));
        v_nl(G.indices(I)) = 1;
        spm_write_vol(h_nl,v_nl);
        if opts.GzipOutput
            gzip(f_nolabel);
            F = appendcleanup(F,f_nolabel);
        end
    else
        f_nolabel = [];
    end
end

%-Cleanup.
for k=1:length(F)
    delete(F{k});
end
end

%==========================================================================
function chk = verifymatch(h1,h2)
chk1 = isequal(h1.dim,h2.dim);
chk2 = all(abs(h1.mat-h2.mat)<1e-6);
chk = and(chk1,chk2);
end

%==========================================================================
function F = appendcleanup(F,f)
if isempty(F)
    F = {f};
else
    F = [
        F
        f
        ];
end
end