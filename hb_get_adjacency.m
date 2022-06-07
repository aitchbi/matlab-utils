function [A,mask] = hb_get_adjacency(mask,conn,varargin)
% HB_GET_ADJACENCY derives an adjacency matrix for a given 2D or 3D mask. 
%
% Inputs:
%   mask: input mask, image or volume, to compute the adjacenct matrix (A)
%   of its pixels/voxels. A will be a square, symmetric matrix, with
%   dimension equal to nnz(mask). mask may alo be absolute address of a
%   nifti file. mask can be binary or grayscale, depending on selected
%   weighting scheme. 
%
%   conn: neighbourhood connectivity for determining pixel/voxel adjacencies. 
%
%   Optional inputs:
%   'weight': edge weighting scheme; default: 'no'.
%
%   'pow': parameter used in weighting scheme: 'yes_power_additive'
% 
%   to be completed...
% 
% Outputs:
%   A: spase adjacency matrix.
%   mask: same as input if dim = length(size(mask)) is 3; if dim=2, then
%   output mask is volume in which two zero images of size input mask have
%   been appended to input mask to make it 3D.
%
% Examples:
%   A=compute_adjacency(mask,26);
%   A=compute_adjacency(mask,6);
%   A=compute_adjacency(mask,26,'weight','no');
%   A=compute_adjacency(mask,26,'weight','yes_power_additive','pow',3);
%
% Acknowledgment:
% Initial version of getedges code by E. Najdenovska & N. Leonardi, 2011.
%
% Hamid Behjat

d = inputParser;
addParameter(d,'weight','no');
addParameter(d,'pow',5);
addParameter(d,'s',[]);
addParameter(d,'distMetric','Eucl');
addParameter(d,'sw',1);
addParameter(d,'sourceFile',[]);
addParameter(d,'verbose',0);
parse(d,varargin{:});
opts = d.Results;

if ischar(mask)
    if opts.verbose
        fprintf(['- source file: ',mask,'\n']); 
    end
    mask = spm_read_vols(spm_vol(mask));
elseif not(isempty(opts.sourceFile))
    if opts.verbose
        fprintf(['- source file: ',opts.sourceFile,'\n']); 
    end
end

if length(size(mask))==2
    fprintf('\n..Input mask was 2D; made 3D.')
    d = zeros([size(mask),3]);
    d(:,:,2) = mask;
    mask = d;
end

[n,m] = getedges(mask,conn);

switch opts.weight
    case 'no'
        if opts.verbose
            fprintf('- weigting: none.\n'); 
        end
        w = ones(length(n),1);
    case 'power_additive'
        % weight defined based on mask values of connected pixels/voxels,
        % as e.g. used in Behjat et al., ISBI-2013 & EMBC-2014.
        pro = 1;
        w = (mask(n).*mask(m)*pro).^opts.pow;
    case 'power_subtractive'
        pro = 1;
        w = (mask(n).*mask(m)*pro).^(-opts.pow);
    case 'nora'
        d = mask(n)-mask(m);
        w = getw(d,opts);
    case 'distance'
        if opts.verbose
            fprintf('- weigting: inverse Euclidean ditance.\n'); 
        end
        [d1,d2,d3] = ind2sub(size(mask),n);
        [d4,d5,d6] = ind2sub(size(mask),m);
        d = sqrt(sum(([d1-d4,d2-d5,d3-d6]).^2,2));
        w = getw(d,opts);
    case 'indexed'
        error('code updated; needs debuging.')
        d = logical(mask(n)-mask(m));
        w = zeros(size(d));
        w(d) = -1; 
end

N = numel(mask);
A = sparse(n,m,w,N,N);
A = A+A';

indices = find(mask);
A = remove_empty_rows_cols(A,indices,'maintain_inds');

end
%==========================================================================
function [X,inds,inds_rmd,c] = remove_empty_rows_cols(X,inds,removeType)

switch removeType
    case 'maintain_inds'
        c=find(~sum(X,1));
        c=c(~ismember(c,inds));
        X(:,c)=[];
        X(c,:)=[];
    case 'update_inds'
        c=find(~sum(X,1));
        X(:,c)=[];
        X(c,:)=[];
        inds_rmd = inds(c);
        inds(c)=[];
end
end
%==========================================================================
function [cInd,nInd] = getedges(mask,conn)

dim = size(mask);
indices = find(mask);
[alli,allj,allk] = ind2sub(dim,indices);

switch conn
    case 6
        nN = 3;
    case 18
        nN = 9;
    case 26
        nN = 13;
    otherwise
        error('Undefined 3D neighbourhood connectivity.');
end

ci = repmat(alli,nN,1);
cj = repmat(allj,nN,1);
ck = repmat(allk,nN,1);

switch conn
    case 6
        ni = [alli  ; alli+1; alli  ];
        nj = [allj+1; allj  ; allj  ];
        nk = [allk  ; allk  ; allk+1];
    case 18
        ni = [alli  ;alli+1;alli+1;alli+1;alli  ;alli  ;alli  ;alli+1;alli+1];
        nj = [allj+1;allj  ;allj-1;allj+1;allj  ;allj+1;allj+1;allj  ;allj  ];
        nk = [allk  ;allk  ;allk  ;allk  ;allk+1;allk-1;allk+1;allk-1;allk+1];
    case 26
        ni = [alli;alli+1;alli+1;alli+1;alli;alli;alli;alli+1;alli+1;alli+1;alli+1;alli+1;alli+1];
        nj = [allj+1;allj;allj-1;allj+1;allj;allj+1;allj+1;allj;allj;allj-1;allj+1;allj+1;allj-1];
        nk = [allk;allk;allk;allk;allk+1;allk-1;allk+1;allk-1;allk+1;allk-1;allk-1;allk+1;allk+1];
end

maskZ = cat(2,mask,zeros(dim(1),1,dim(3)));
maskZ = cat(1,maskZ,zeros(1,dim(2)+1,dim(3)));
maskZ = cat(3,maskZ,zeros(dim(1)+1,dim(2)+1,1));

valid = (ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2) & nk>=1 & nk<=dim(3));
ni = ni(valid);
nj = nj(valid);
nk = nk(valid);
ci = ci(valid);
cj = cj(valid);
ck = ck(valid);

tt = sub2ind(size(maskZ),ni,nj,nk);
ee = maskZ(tt);
valid = logical(ee);
ni = ni(valid);
nj = nj(valid);
nk = nk(valid);
ci = ci(valid);
cj = cj(valid);
ck = ck(valid);

cInd = sub2ind(dim,ci,cj,ck);
nInd = sub2ind(dim,ni,nj,nk);
end
%==========================================================================
function w = getw(d,opts)
if isempty(opts.s)
    opts.s = mean(abs(d))*opts.sw;
end
switch opts.distMetric
    case 'Gaus'
        w = exp(-d.^2/(2*opts.s^2));  % Gaussian similarity
    case 'adjEucl'
        w = 1./(1+d.^2/(2*opts.s^2)); % adjusted Euclidean similarity
    case 'Eucl'
        w = 1./d;                     % Euclidean similarity
end
end
