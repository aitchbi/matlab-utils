function X = hb_zscore(X,n,I)
%HB_ZSCORE z-scores input matrix X; either the rows or columns.
%
% Inputs:
%   X: a matrix.
%   n: n=1: z-score rows, n=2: z-score columns.
%   I: logical vector of length size(X,2), specifying which nodes to
%   operate on and which ones to skip. 
%
% Hamid Behjat 

if ~exist('I','var') || isempty(I)
    I = true(1,size(X,n));
end

assert(length(size(X))==2);
assert(any(n==[1 2]));
if n==1
    X = X';
end

% skip some columns
Xskip = X(:,not(I));
X = X(:,I);

% z-score
N = size(X,2); 
L = 1e5; 
K = ceil(N/L); 
for k=1:K
    if k~=K
        kk = (k-1)*L+1:k*L;
    else
        kk = (k-1)*L+1:N;
    end
    X(:,kk) = zscore(X(:,kk));
end

% put back skiped columns
if not(isempty(Xskip))
    d = zeros(size(X,1),length(I));
    d(:,not(I)) = Xskip;
    d(:,I) = X;
    X = d;
end

if n==1
    X = X';
end
end