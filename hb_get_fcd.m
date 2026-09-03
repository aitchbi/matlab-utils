function FCD = hb_get_fcd(X, W)
% HB_GET_FCD derives a functional connectivity dynamics (FCD) matrix based
% on input functional time series using a sliding-window approach.
%
%   FCD = HB_GET_FCD(X, W)
%
% inputs
%   X : NxT matrix of functional (e.g., fMRI) time series where N is the
%       number of regions and T is the number of time points.
%   W : Window length in time points (e.g., use W = 83 for HCP data, which
%   is equivalen to ~60 seconds).
%
% outputs
%   FCD : MxM matrix of window-to-window similarity of functional
%         connectivity patterns, where M = T - W + 1.
%
% method
%   as in https://doi.org/10.1073/pnas.2318641121 
%   • slide a window of length W across the T samples producing M windows.
%   • for each window, compute the region-by-region FC matrix (Pearson r).
%   • vectorize the upper triangle (excluding the diagonal) of each FC.
%   • correlate these vectors across all windows to obtain the MxM FCD.
%
% example
%   % given an X of size N x 1200 and W = 83, FCD will be a 1118 x 1118
%   matrix.
%
%   FCD = hb_get_fcd(X, 83);
%
% HB

%-check input.
%--------------------------------------------------------------------------
if nargin < 2
    error('hb_get_fcd:NotEnoughInputs', 'Provide X (NxT) and W (window length).');
end

if ~isnumeric(X) || ndims(X) ~= 2
    error('hb_get_fcd:BadX', 'X must be a numeric 2D matrix of size [N x T].');
end

[N, T] = size(X);

if ~isscalar(W) || ~isfinite(W) || W ~= floor(W) || W < 2 || W > T
    error('hb_get_fcd:BadW', 'W must be an integer in [2, T].');
end

%-prepare indices and storage.
%--------------------------------------------------------------------------
M = T - W + 1; % number of sliding windows
L = N * (N - 1) / 2; % number of unique FC edges (upper triangle, no diag)
triu_idx = find(triu(true(N), 1));

V = zeros(L, M, 'like', X); % columns are vectorized FCs per window

%-compute FC for each window and vectorize.
%--------------------------------------------------------------------------
for m = 1:M
    t0 = m;        
    t1 = m + W - 1;
    Xw = X(:, t0:t1); % N x W
    
    % corrcoef input format: 
    % variable/regions-in-columns, observations/timpepoints-in-rows
    R = corrcoef(Xw'); % N x N

    % upper triangle, excluding diagonal
    V(:, m) = R(triu_idx); % L x L
end

%-correlate FC vectors across windows to form FCD.
%--------------------------------------------------------------------------
% corrcoef input format: 
% variables/windows-in-colum, observations/edges-in-rows
FCD = corrcoef(V); % M x M

%-ensure numerical symmetry & bounds.
%--------------------------------------------------------------------------
FCD = max(min((FCD + FCD.') / 2, 1), -1); % handles tiny asymmetries
end