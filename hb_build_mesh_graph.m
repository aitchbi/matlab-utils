function [A, I, mask] = hb_build_mesh_graph(mask, conn)
% HB_GET_MESH_GRAPH builds adjacency matrix of a mesh graph from a 3D mask.
% 
% inputs: 
%   mask: binary numeric 3D array or a nifti file; if the latter, file will
%   be made binary if after loading if not already so.
%
%   conn: neighbourhood connectivity in 3D; 6, 8, or 26.
%
% output:
%   A: mesh graph adjacency matrix; a sparse matrix.  
%
% h behjat

if ischar(mask)
    d1 = endsWith(mask, '.nii');
    d2 = endsWith(mask, '.nii.gz');
    assert(d1||d2);
    if d1
        CleanUp = 0;
    else
        gunzip(mask);
        CleanUp = 1;
    end
    mask = spm_read_vols(spm_vol(mask));
    if CleanUp
        delete(strrep(mask, '.nii.gz', '.nii'));
    end
end

dim = size(mask);

I = find(mask);

[ii, jj, kk] = ind2sub(dim, I);

switch conn
    case 6
        nN = 3;
    case 18
        nN = 9;
    case 26
        nN = 13;
end

cx = repmat(ii, nN, 1);
cy = repmat(jj, nN, 1);
cz = repmat(kk, nN, 1);

switch conn
    case 6
        nx = [ii  ; ii+1; ii  ];
        ny = [jj+1; jj  ; jj  ];
        nz = [kk  ; kk  ; kk+1];
    case 18
        nx = [ii  ;ii+1;ii+1;ii+1;ii  ;ii  ;ii  ;ii+1;ii+1];
        ny = [jj+1;jj  ;jj-1;jj+1;jj  ;jj+1;jj+1;jj  ;jj  ];
        nz = [kk  ;kk  ;kk  ;kk  ;kk+1;kk-1;kk+1;kk-1;kk+1];
    case 26
        nx = [ii;ii+1;ii+1;ii+1;ii;ii;ii;ii+1;ii+1;ii+1;ii+1;ii+1;ii+1];
        ny = [jj+1;jj;jj-1;jj+1;jj;jj+1;jj+1;jj;jj;jj-1;jj+1;jj+1;jj-1];
        nz = [kk;kk;kk;kk;kk+1;kk-1;kk+1;kk-1;kk+1;kk-1;kk-1;kk+1;kk+1];
end

maskZ = cat(2, mask, zeros(dim(1), 1, dim(3)));
maskZ = cat(1, maskZ, zeros(1, dim(2)+1, dim(3)));
maskZ = cat(3, maskZ, zeros(dim(1)+1, dim(2)+1, 1));

valid = (nx>=1 & nx<=dim(1) & ny>=1 & ny<=dim(2) & nz>=1 & nz<=dim(3));
nx = nx(valid);
ny = ny(valid);
nz = nz(valid);
cx = cx(valid);
cy = cy(valid);
cz = cz(valid);

tt = sub2ind(size(maskZ),nx,ny,nz);
ee = maskZ(tt);
valid = logical(ee);
nx = nx(valid);
ny = ny(valid);
nz = nz(valid);
cx = cx(valid);
cy = cy(valid);
cz = cz(valid);

cInd = sub2ind(dim, cx, cy, cz);
nInd = sub2ind(dim, nx, ny, nz);

% build sparse matrix
N = prod(dim);
A = sparse([cInd, nInd], [nInd, cInd], ones(1, 2*numel(nx)), N, N);

% remove empty rows/columns
c = find(~sum(A, 1));
c = c(~ismember(c, I));
A(:,c) = [];
A(c,:) = [];
end


