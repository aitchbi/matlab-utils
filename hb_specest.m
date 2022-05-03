function [counts,bcents,sts,bedges,bwidth,funcPath] = hb_specest(L,R,N,todo,th,sp,tt)
%HB_SPECEST estimates the distribution of eigenvalues of Laplacian matrix 
% L at N bins within range R, where the bin widths are fixed and equal to 
% (R(2)-R(1))/N. If R(2) > maximum eigenvalue of L (lmax), counts for bins 
% that fall above lmax are set to nan.
%
% Inputs:
%   L   : graph Laplacian matrix. 
%   R   : 1x2 scalar array, spectral range.
%   N   : number of subbands.
%   todo: (optional) 1xN logical array specifying which of the N subbands 
%           to compute; [default: true(1,N)] 
%   th  : (optional) scalar in interval [0, 0.5]; this is the pivot 
%           tolerance used for ldl factorization in ldl.m; [default: 0.01]   
%   sp  : (optional) show progress; [default: true]
%   tt  : (optional) tic-toc. 
%               
% Output:
%   counts : 1xN array of eigenvalue counts.
%   bcents : 1xN array of bin centers.
%   sts    : estimation status.
%   bedges : 1x(N+1) array of bin edges; bedges(1)=R(1), bedges(end)= R(2).
%   bwidth : bin width.
%   funcPath: path to this function. 
%
% Examples:
% C = hb_specest(L,[0,2],100); % distribution of eigenvals upto 2, 100 bins
% C = hb_specest(L,[0,2],100,[true,false(1,N-1)]); % compute only 1st bin
% C = hb_specest(L,[0,.1],5); % distrib. of eigenvals upto 0.1, 5 bins
% C = hb_specest(L,[0,.1],1); % total # of eigenvalues below 0.1
% C = hb_specest(L,[.1,1],5); % distrib. of eigenvals within [.1,1], 5 bins
% C = hb_specest(L,[.1,1],1); % total # of eigenvals within [.1,1]
% C = hb_specest(L,[0,2],100,[],0.01); % modify pivot tolerance in ldl.m 
% C = hb_specest(L,[0,2],200,[],[],false); % display off for batch runs. 
% [C,~,sts] = hb_specest(L,[0,2],200); % sts: estimate ok?

% Hamid Behjat

if ~exist('todo','var')||isempty(todo)
    todo = true(1,N);
else
    assert(length(todo)==N);
end
if ~exist('th','var')||isempty(th)
    th = .01;
end
if ~exist('sp','var')||isempty(sp)
    sp = true; 
end
if ~exist('tt','var')||isempty(tt)
    tt = false; 
end

lmax = eigs(L,1,'largestreal');

if R(1)>=lmax, error(''); end
if R(1)<0, error(''); end
if R(1)>=R(2), error(''); end

bwidth = (R(2)-R(1))/N; % bin width
bedges = R(1):bwidth:R(2); % bin edges
bcents = bedges(1:end-1)+bwidth/2;
Ne = length(bedges);

if R(2)>lmax
    iu = find(bedges==lmax);
    if isempty(iu)
        iu = find(bedges>lmax,1,'first');
        bedges(iu) = lmax;
    end
else
    iu = length(bedges);
end
if iu<=N
    todo(iu:end) = false;
end

I = find(todo); % lower bin edges
Ni = length(I);

if sp
    fprintf('\n..Estimating spectrum..');
    fprintf('\n  # of bins: %d',Ni);
    fprintf('\n  bin width: %d',bwidth);
end

counts = nan(Ne,1);
if bedges(1)==0
    counts(1) = 0;
end
if bedges(end)==lmax
    counts(Ne) = size(L,1);
end

for i = 1:Ni
    if tt, tic; end
    k = I(i); % lower bin edge
    if isnan(counts(k))
        counts(k) = ldlfact(L,bedges(k),th);
    end
    k = I(i)+1; % upper bin edge
    if isnan(counts(k))
        counts(k) = ldlfact(L,bedges(k),th);
    end
    if sp, progress(i,Ni); end
    if tt, toc, end
end
counts = diff(counts(:))';

if any(counts<0)
    sts = ['wrong estimate in ',num2str(nnz(counts<0)),' of the bins.'];
else
    sts = 'ok estimate.';
end
funcPath = mfilename('fullpath');
end
%==========================================================================
function c = ldlfact(L,be,th)
P = symamd(L);
T = L-be*speye(size(L,1));
[~,D,~] = ldl(T(P,P),th);
c = sum(diag(D)<0);
end
%==========================================================================
function progress(n,m)
tag = '  ';
l = length(num2str(m));
if n==1
    fprintf(['\n',tag]);
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,m)'])
end