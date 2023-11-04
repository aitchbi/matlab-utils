function [v_o, h_o, C, cost, lbls, BgType] = hb_nii_match_labels(f_i, f_r, f_o, varargin)
% HB_NII_MATCH_LABELS relables segmentation file f_i to have labels that
% optimaly match those in f_r. By default, dice score is used as the
% similarity measure for pairs of segments and the Hungarian matching
% algorith (munkres) is used to perform assignments.
% 
%-Exp 1: pecify output file name
% hb_nii_match_labels(f_i, f_r, 'Path/to/f_o.nii');
%
%-Exp 2: overwrite input file
% hb_nii_match_labels(f_i, f_r, f_i); 
%
% Dependencies:
%   hb_nii_load.m
%   hb_nii_write.m
%   SPM12 package
%   munkres.m
%
% Hamid Behjat

d = inputParser;
addParameter(d,'MatchingMethod', 'munkres'); 
addParameter(d,'SimilarityMeasure', 'dice'); % 'dice', 'centroid'
addParameter(d,'BackgroundVoxelLabel', []);
addParameter(d,'DontWrite', false);
parse(d,varargin{:});
opts = d.Results;

%-[a] Handle inputs.
[v_i, v_r, h_o] = handleinputs(f_i, f_r, f_o);

%-[b] Handle NaNs & ensure non-negative labels.
[v_i, BgType] = fixbg(v_i);
[v_r, d]      = fixbg(v_r);
assert(isequal(BgType, d));

%-[c] Ensure matching label sets.
lbls = unique(v_i(:));
d = unique(v_r(:));
assert(isequal(d, lbls), 'Labels of input & ref must match.');
switch BgType
    case '0'
        assert(lbls(1)==0); 
        lbls = lbls(2:end); % skip background (0).
    case 'NaN'
        assert(lbls(1)==-1); 
        lbls = lbls(2:end); % skip background (-1).
        % lbls(2) can be 0 or above, since 0 is not background. 
end
N = length(lbls);

%-Biuld cost matrix.
C = zeros(N);
% C(i,j): cost of assigning lbls(j) in v_i lbls(i) in v_r

switch opts.SimilarityMeasure
    case 'centroid'

        % parcel centroids 
        xyz_i = get_centroids(v_i, lbls);
        xyz_r = get_centroids(v_r, lbls);

        % Euclidean distance between centroids
        for k=1:N
            d = repmat(xyz_r(k,:),N,1)-xyz_i;
            d = d.^2;
            d = sqrt(sum(d,2));
            C(k,:) = d';
        end

    case 'dice'
        S = size(v_i);
        v0 = zeros(S);
        for a=1:N
            d1 = v0;
            d1(v_r==lbls(a)) = 1;
            for b=1:N
                d2 = v0;
                d2(v_i==lbls(b)) = 1;
                C(a,b) = nnz(d1 & d2) / nnz(d1 | d2);
            end
        end
        C = 1-C;
end

%-Do matching.
switch opts.MatchingMethod
    case 'munkres'
        [assignment, cost] = munkres(C);
    otherwise

end

%-Build re-labeled file.
switch BgType
    case '0'
        v_o = zeros(size(v_i));
    case 'NaN'
        v_o = nan(size(v_i));
end
for k=1:N
    I = find(v_i==lbls(assignment(k)));
    v_o(I) = lbls(k);
end

%-Write relabelled file.
if opts.DontWrite
    % just returning volume; no writing. 
else
    hb_nii_write(h_o, v_o);
end
end

%==========================================================================
function [v_i, v_r, h_o] = handleinputs(f_i, f_r, f_o)
[v_i, h_i] = hb_nii_load(f_i);
[v_r, h_r] = hb_nii_load(f_r);
assert(all(isequal(h_i.dim,h_r.dim)));
assert(all(all(abs(h_i.mat-h_r.mat)<1e-6)));
h_o = h_i;
h_o.fname = f_o;
end

%==========================================================================
function [v, BgType] = fixbg(v)
p = 30; % percent
n = numel(v)*p/100;

d = nnz(isnan(v));
if d==0
    d = nnz(v==0);
    assert(d > n);
    BgType = '0';
else
    if d>n
        % more than p% of volume is NaN, thus, background is NaN
        assert(not(any(v(:)<0))); % ensure all labels are positive
        v(isnan(v)) = -1;         % set background to -1.
        BgType = 'NaN';
    else
        msg = 'Fishy NaNs in file that are seemingly not the background.';
        error(msg);
    end
end
end

%==========================================================================
function xyz = get_centroids(v, lbls)
S = size(v);
N = length(lbls);
D = length(S);
xyz = zeros(N, D);
for k=1:N
    d = zeros(S);
    l = lbls(k);
    d(v==l) = 1;
    switch D
        case 2
            d = regionprops(d, 'Centroid');
        case 3
            d = regionprops3(d, 'Centroid');
    end
    xyz(k,:) = d.Centroid;
end
end