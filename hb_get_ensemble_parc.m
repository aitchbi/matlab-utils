function [v_pc, v_pb, labels, BgLabel] = hb_get_ensemble_parc(F, varargin)
% HB_GET_ENSEMBLE_PARC generates an ensemble parcellation and probability
% map based on an ensemble set of parcellations.
%
% Inputs:
%   F: a cell array of parcellation (segmentation) filenames from which a
%   single parcellation map is to be extracted. By default, F{1} is treated
%   as the reference file, to which the other files are initially
%   label-matched via the Hungarian matching algorithm; another file can be
%   specified as the reference via optional Name-Value input argument
%   'WhichFileIsReference'. The label matched files are then aggregated to
%   obtain a single parcellation file (v_pc) in which each voxel is
%   assigned the label which is most seen across files, at least 50% of the
%   files. In doing so, some voxels, manily at boundaries of parcels may
%   not be assigned a label. The ensemble parcellation is accompanied with
%   a probability map (v_pb), showing the probability of the assignment of
%   each voxel; a probablity of 1 means that a voxel was label the same
%   across all files, whereas 0.2 means on 20% of the files.
%
% Note: if value of background voxels not specifed via the optional
% Name-Value input argument 'BackgroundVoxelLabel', the alogorithm will
% determine it, considering that it is either 0 or NaN. Furthermore, the
% assumption is that the file does not contain NaNs unless as background.
%
% Dependencies:
%   hb_nii_load.m
%   hb_nii_write.m
%   hb_nii_match_labels.m
%   SPM12 package
%   munkres.m
%
% Hamid Behjat

d = inputParser;
addParameter(d,'MatchingMethod', 'munkres');
addParameter(d,'SimilarityMeasure', 'dice'); % 'dice', 'centroid'
addParameter(d,'UnassignedVoxelLabel', []);  % an integer or NaN
addParameter(d,'BackgroundVoxelLabel', []);  % an integer or NaN
addParameter(d,'WhichFileIsReference', []);  % an integer in [1 length(F)] 
addParameter(d,'OutputFilenameParc', []);    % pathname of file to write
addParameter(d,'OutputFilenameProb', []);    % pathname of file to write

parse(d,varargin{:});
opts = d.Results;

%-Reference file.
if isempty(opts.WhichFileIsReference)
    f_r = F{1};
else
    f_r = F{opts.WhichFileIsReference};
end
[v_r, h_r] = hb_nii_load(f_r);

%-Get label indices.
I_nan = find(isnan(v_r));
if isempty(opts.BackgroundVoxelLabel)
    if isempty(I_nan)
        BgLabel = '0';
    else
        BgLabel = 'NaN';
    end
else
    BgLabel = opts.BackgroundVoxelLabel;
end
d = find(v_r);
I_lbl = setdiff(d, I_nan); % non-zero/NaN
if strcmp(BgLabel, 'NaN')
    % 0 may be a label
    d = find(v_r==0);
    I_lbl = union(I_lbl, d);
end

%-Match labels.
N_fls = length(F);
N_vox = length(I_lbl);
L = zeros(N_fls, N_vox);
for k=2:N_fls
    f = F{k};
    [d, ~, ~, ~, lb, bg] = ...
        hb_nii_match_labels(...
        f, f_r, [],...
        'SimilarityMeasure', opts.SimilarityMeasure, ...
        'MatchingMethod', opts.MatchingMethod,...
        'BackgroundVoxelLabel', BgLabel,...
        'DontWrite', true);
    assert(isequal(bg, BgLabel));
    L(k,:) = d(I_lbl(:))';
    if k==2
        labels = lb;
    else
        assert(isequal(lb, labels));
    end
end
L(1,:) = v_r(I_lbl(:))';

%-Find modes.
[l, m, c] = mode(L); % NOTE 1

%-Build ensemble.
d = size(v_r);
switch BgLabel
    case '0'
        v_pc = zeros(d);
        v_pb = zeros(d);
    case 'NaN'
        v_pc = nan(d);
        v_pb = nan(d);
end
v_pc(I_lbl) = l;
v_pb(I_lbl) = m/N_fls;

%-Handle unassigned voxels.
d = cellfun(@length, c);
I_ua = I_lbl(d>1); % NOTE 2
if not(isempty(I_ua))
    if isempty(opts.UnassignedVoxelLabel)
        switch BgLabel
            case '0'
                l_ua = 0;
            case 'NaN'
                l_ua = NaN;
        end
    else
        l_ua = opts.UnassignedVoxelLabel;
        d1 = isinteger(l_ua);
        d2 = isnan(l_ua);
        d = 'Integer or NaN expected as label.';
        assert(or(d1, d2), d);
    end
    v_pc(I_ua) = l_ua; % NOTE 3
end

%-Write files.
if not(isempty(opts.OutputFilenameParc))
    % ensemble parcellation
    h = struct;
    h.dt    = h_r.dt;
    h.dim   = h_r.dim;
    h.mat   = h_r.mat;
    h.fname = opts.OutputFilenameParc;
    hb_nii_write(h, v_pc);
end

if not(isempty(opts.OutputFilenameProb))
    % ensemble probability map
    h = struct;
    h.dt    = [64 0];
    h.dim   = h_r.dim;
    h.mat   = h_r.mat;
    h.fname = opts.OutputFilenameProb;
    hb_nii_write(h, v_pb);
end

end

%-NOTES.
%
% NOTE 1
% l: most frequent label per voxl
% m: frequency of most repeated label
% c: all labels that have same frequency as l
%
% NOTE 2
% Ambigious voxels having two or more equaly frequent labels. These voxels
% are not assigned a parcellation label. Instead, by default, these voxels
% are assigned the value of the background voxels, or a specific label can
% be provide by the optional Name-Value argument 'UnassignedVoxelLabel'.
%
% NOTE 3
% Only within v_prc we assign a specific label to unassigned voxels. Since
% v_prb is a probability map, we retain the probability even at unassigned
% voxels.
