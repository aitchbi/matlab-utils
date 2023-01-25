function X = hb_detrend(X,n,polyords,I)
%HB_DETREND detrends each row (n=1) or column (n=2) of matrix X, chunk by
% chunk if too many rows or columns, respectively. 
%
% Inputs:
%   X: matrix.
%   n: 1 (row) or 2 (column).
%   polyords: a non negative integer scalar or vector; polynomial orders.
%   I: inclusion list; logical array, specifying which row/columns to
%   do/skip.
%
% Outputs: 
%   X: detrended matrix.
%
% Examples:
%   X = hb_detrend(X,1,1);     % each row: remove linear trend.
%   X = hb_detrend(X,1,[0 1]); % each row: remove mean & linear trend.
%   X = hb_detrend(X,1,[1 2]); % each row: remove linear & quad trends.
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

% skip some data
Xskip = X(:,not(I));
X = X(:,I);

N = size(X,2);
L = 1e5;
K = ceil(N/L);
for k=1:K
    if k~=K
        kk = (k-1)*L+1:k*L;
    else
        kk = (k-1)*L+1:N;
    end
    for iOrd=1:length(polyords)
        ord = polyords(iOrd);
        X(:,kk) = detrend(X(:,kk),ord);
    end
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