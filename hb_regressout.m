function X = hb_regressout(X,R,runpar,Np,I)
%HB_REGRESSOUT Regress out effect from regressors given as columns in
%matrix R from each column in matrix X. 
%
% Inputs:
%   X: T*N matrix, where N is the number of vectors to process.
%   R: T*M matrix, where M is the number of regressors.
%   runpar: (optional) run in parallel over columns in X. (default: false)
%   Np: (optional) number of parallel workers to use; default: use max number
%   available as determined by gcp. 
%   I: logical vector of length size(X,2), specifying which nodes to
%   operate on and which ones to skip. 
%
% Outputs:
%   X: X after regressing out effects from regressors.  
%
% Hamid Behjat

if ~exist('runpar','var') || isempty(runpar)
    runpar = false; % spmd didn't work so well with HCP rest as X
end
if ~exist('Np','var')
    Np = 2; % be cautios of large memory usage for large X 
end
if ~exist('I','var') || isempty(I)
    I = true(1,size(X,2));
end

assert(length(size(X))==2);
assert(size(X,1)==size(R,1),'number of rows in X and R shoud be equal.');

% leave out some columns
Xskip = X(:,not(I));
X = X(:,I);

% regress out
N = size(X,2);
if ~runpar
    for k=1:N
        b = regress(X(:,k),R);
        r = sum(repmat(b(:)',size(X,1),1).*R,2);
        X(:,k) = X(:,k)-r;
    end
    %X = dorun(X,R); computationally more efficient by not calling function
else
    % Split across workers.
    [x,p] = hb_spmd_prepare([{[]},{R},{1:N}],[0,0,1],Np);
    
    % Split X to make run memory-safe.
    x{1} = cell(1,p.NumWorkers);
    for k=1:p.NumWorkers
        x{1}{k} = X(:,x{3}{k});
    end
    
    % Run on workers. 
    spmd
        y = dorun(x{1}{labindex},x{2});
    end
    
    % Combine results from workers.
    for k=1:p.NumWorkers
        X(:,x{3}{k}) = y{k};
    end
end

% put back skipped columns
if not(isempty(Xskip))
    d = zeros(size(X,1),length(I));
    d(:,not(I)) = Xskip;
    d(:,I) = X;
    X = d;
end
end

function X = dorun(X,R)
T = size(X,1);
N = size(X,2);
for k=1:N
    b = regress(X(:,k),R);
    r = sum(repmat(b(:)',T,1).*R,2); 
    X(:,k) = X(:,k)-r;
end
end