function [w, e, w_ns, Y] = hb_get_fine_ee_warping(A,B,varargin)

d = inputParser;
addParameter(d,'SmoothSpan', []);
addParameter(d,'TypeB', 'GFT'); % 'GFT' or 'energy'
parse(d,varargin{:});
opts = d.Results;

assert(isvector(A));
assert(isvector(B));

K = length(A);
assert(K==length(B));

switch opts.TypeB
    case 'GFT'
        % convert GFT coefficients to energy
        B = cellfun(@(x) x.^2, B, 'UniformOutput', false); 
    case 'energy'
        % B already is energy 
end

lmax = min(cellfun(@max, A));

M = round(mean(cellfun(@length, A)));

e = linspace(0,lmax,M);

BC = cellfun(@(x) cumsum(x,1), B, 'UniformOutput', false);

Y_sum = zeros(M,1);
Y     = cell(K,1);
for k=1:K
    J = size(B{k},2);
    Y{k} = zeros(M,J);
    Ak = A{k}(:);
    for j=1:J
        Y{k}(:,j) = interp1(Ak, BC{k}(:,j), e);
    end
    if k==1
        J_sum = J;
        Y_sum = sum(Y{k},2);
    else
        J_sum = J_sum + J;
        Y_sum = Y_sum + sum(Y{k},2);
    end
end
w_ns = Y_sum/J_sum;
w_ns = w_ns./w_ns(end);

%-Smooth.
L = opts.SmoothSpan;
if isempty(L)
    L = lmax/0.05;
elseif L<1
    % L: a spectral length
    d = e(2)-e(1);
    L = round(L/d);
else
    assert(isinteger(L));
end
p = 0.1; % p% of input spectral range
n = p*M;
assert(L<n);
w = smooth(w_ns, L);
w = w./w(end);
end