function [y,L,r,o] = hb_gsig_filt(A,x,g,varargin)
% HB_GSIG_FILT filters a given graph signal using a given
% set of spectral graph kernels using Chebyshev polynomial approximation,
% i.e., vertex-domain filtering without diagonalizing the shift operator
% (e.g. graph Laplacian matrix). The current implementation is only
% applicable to the following shift operators:
%
%   .Combinatorial Laplacian
%   .Symmetric Normalized Laplacian [default] 
%   .Random-walk Laplacian 
%
% Inputs:
%   A: graph adjacency matrix.
%   x: graph signals, one per column. 
%   g: cell array of spectral kernels. 
%
% Outputs: 
%   y: filtered graph signals, cell array of length size(x,2), each
%   containing a size(x,1) x length(g) matrix.
%   r: spectrum lower and upper range. 
%   o: see above, if not given as input. 
%
% Examples:
% y = hb_gsig_filt(A,x,g);
% y = hb_gsig_filt(___,'ChebypolOrds',o);
% y = hb_gsig_filt(___,'ChebypolOrds',o,'KernelsShouldFormTightFrame',1);
%
% Dependencies:
% .SGWT: https://wiki.epfl.ch/sgwt/documents/sgwt_toolbox-1.02.zip
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_cheby_order_est.m
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_laplacian.m
%
% Hamid Behjat

%-Process inputs. 
opts = processinputs(A,x,g,varargin,inputParser);

%-Laplacian etc. 
[L, r] = getlaplacian(A, opts);

%-Verify/estimate Chebyshev polynomial order.
o = getchebyords(opts);

%-Compute Chebyshev polynomial coefficients.
c = getchebucoeffs(g, o, r);

%-Filter graph signals.
y = dofilt(x,L,c,r,opts);

end

%==========================================================================
function opts = processinputs(A,x,g,varinputs,p)

% verify required inputs 
assert(size(x,1)==size(A,1));
assert(iscell(g));

% handle optional inputs
f1 = @(x) assert(ischar(x) && ismember(x,{'cb','sn','rw'}));
f2 = @(x) assert(isnumeric(x) && length(x)==length(g));
f3 = @(x) assert(all([isnumeric(x),all(x>0),length(x)==2]));
f4 = @(x) assert(all([isnumeric(x),x>0,length(x)==1]));
f5  = @(x) assert(islogical(x) || ismember(x,[0 1]));

%p = inputParser;
addParameter(p,'ShifOperator', 'sn', f1);
addParameter(p,'ChebypolOrds', [], f2); % [1]
addParameter(p,'ChebypolOrdSearchRange', [10 500], f3); % [2]
addParameter(p,'ChebypolOrdErrorTolKernel', 1e-4, f4); % [3]
addParameter(p,'ChebypolOrdErrorTolTightframe', 1e-4, f4); % [4]
addParameter(p,'KernelsShouldFormTightFrame', false, f5); % [5]
addParameter(p,'Verbose', false, f5); 
parse(p,varinputs{:});
opts = p.Results;
% [2]: only applicable if [1] is empty
% [3]: only applicable if [1] is empty
% [4]: only applicable if [1] is empty
% [5]: only applicable if [1] is empty & [5] is true. 
%
% [1]: structure with fields 'kernel' and 'tightframe' each of length(g),
% specifying the Chebyshev polynomial order to use for appriximating each
% kernel in g. If [], polynomial orders will be estimated.
end

%==========================================================================
function [L, arange] = getlaplacian(A, opts)
so = opts.ShifOperator;
switch so
    case 'sn'
        try 
            L = hb_laplacian(A, so);
        catch 
            L = sgwt_laplacian(A, 'opt', 'normalized');
        end
    otherwise
        L = hb_laplacian(A, opts.ShifOperator);
end
switch so
    case {'cb', 'sn', 'rw'}
        lmax = getlmax(L);
        switch so
            case 'cb'
                % no limit on upper bound
            case {'sn', 'rw'}
                % upper bound of lmax is 2
                lmax = min(lmax, 2);
        end
        arange = [0 lmax];
    otherwise
        error('extend');
end
end

%==========================================================================
function cOrds = getchebyords(opts)
cOrds = opts.ChebypolOrds;
if isempty(cOrds)
    error('extend');
    if 1
        opts_cOrds = [];
    else
        opts_cOrds = struct;
        opts.ChebypolOrdSearchRange;
    end
    [cOrds,e,G,Gc] = hb_cheby_order_est(g,arange,opts_cOrds);
end
end

%==========================================================================
function c = getchebucoeffs(g, cOrds, arange)
J = length(g);
c = cell(1,J);
for k = 1:J
    d = cOrds(k);
    c{k} = sgwt_cheby_coeff(g{k},d,d+1,arange);
end
end

%==========================================================================
function y = dofilt(x,L,c,r,opts)
if opts.Verbose
    fprintf('\n.Filtering signals.. [%d signals, %d filters]',...
        size(x,2),...
        length(c));
end
y = sgwt_cheby_op(x, L, c, r);
end

%==========================================================================
function lmax = getlmax(L)
opts = struct('tol', 5e-3, 'p', 10, 'disp', 0);
lmax = eigs(L, 1, 'lm', opts);
lmax = lmax*1.01;
end
