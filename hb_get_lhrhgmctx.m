function [f_o, G] = hb_get_lhrhgmctx(G, varargin)
%
% Inputs:
%   G: voxel-wise graph structure.
%
% Name-Value Pair Arguments:
%   DirSubjects: root forlder where your subjects are saved.
%
%   GraphRefNii: a reference nifti in the space where the graph was build;
%   if not given, will use G.f.mask if availble. 
%
% Outputs:
%   f_o: gray matter cerebral cortex mask associated to the graph, wherein
%   left and right hemispheres are marked. In particular, binary mask file
%   G.f.mask_gm_ctx is loaded and then split into left (label: 3) and right
%   hemisphere (label:42); fishy voxels are are given label=50. If
%   G.f.mask_gm_ctx does not exist due to incorrect file path*,
%   <DirSubjects> will be used to correct the filepath; *e.g. graph
%   originally generated on another machine and thus filepaths not matching
%   present machine paths.
%
%   G: Updated graph structure, with fields specifying indices of left and
%   right hemispheres, and unassigned ones.
%
% Dependencies:
% .SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_nii_reslice.m 
% .https://github.com/aitchbi/matlab-utils/blob/main/spm_modified
%
% HB

d = inputParser;
addParameter(d,'DirSubjects', []);
addParameter(d,'OutputFile',[]);
addParameter(d,'GraphRefNii',[]); % if set, also should set OutputFile
addParameter(d,'OverwriteExistingFile',false);
parse(d,varargin{:});
opts = d.Results;

if isempty(opts.GraphRefNii)
    f_gmctx = G.f.mask_gm_ctx; % GM ctx mask
    TmpGmCtx = false;
    N_gmctx = length(G.indices_gm_ctx);
    if isempty(opts.OutputFile)
        f_o = strrep(f_gmctx, '.nii', '_lhrh.nii');
    else
        f_o = opts.OutputFile;
    end
    [OutputFileExists, f_o] = outputfileexists(f_o, opts);
else
    assert(not(isempty(opts.OutputFile)));
    f_o = opts.OutputFile;
    [OutputFileExists, f_o] = outputfileexists(f_o, opts);
    if not(OutputFileExists)
        f_ref = opts.GraphRefNii;
        if exist(f_ref, 'file')
            f_refgz = [];
        elseif exist([f_ref, '.gz'], 'file')
            f_refgz = [f_ref, '.gz'];
            gunzip(f_refgz);
        end
        h_ref = spm_vol(f_ref);
        if not(isempty(f_refgz)) % cleanup gunziped file
            delete(f_ref);
        end
        if length(h_ref)>1
            h_ref = h_ref(1);
        end
        f_gmctx = strrep(f_o, '.nii', '_tmp_gmctx.nii');
        TmpGmCtx = true;
        h_gmctx = struct;
        h_gmctx.fname = f_gmctx;
        h_gmctx.dim = h_ref.dim;
        h_gmctx.mat = h_ref.mat;
        h_gmctx.dt = [16 0];
        v_gmctx = zeros(h_gmctx.dim);
        v_gmctx(G.indices_gm_ctx) = 1;
        spm_write_vol(h_gmctx,v_gmctx);
        N_gmctx = length(G.indices_gm_ctx);
    end
end

if OutputFileExists
    v_gmctxlhrh = spm_read_vols(spm_vol(f_o));
    UniqueLabels = sort(unique(v_gmctxlhrh(:)));
else
    % GM ctx mask
    [h_gmctx, v_gmctx] = getgmctx(G, f_gmctx, opts.DirSubjects);

    % reference ribbon
    v_rib = getrib(G, f_gmctx, opts.DirSubjects);

    % lh/rh labeled GM ctx mask
    v_gmctxlhrh = getlhrhvol(v_gmctx, v_rib, N_gmctx);
    UniqueLabels = writelhrhvol(v_gmctxlhrh, v_gmctx, N_gmctx, h_gmctx, f_o);

    % cleanup
    if TmpGmCtx
        delete(f_gmctx);
    end
end

% update G
G = updateG(G, v_gmctxlhrh, UniqueLabels);
end

%==========================================================================
function [FileExists, f] = outputfileexists(f, opts)
fgz = [f, '.gz'];

if not(opts.OverwriteExistingFile)
    if exist(f, 'file')
        fprintf('\n.File exists: %s\n', f);
        FileExists = true;
    elseif exist(fgz, 'file')
        fprintf('\n.Gziped file exists: %s\n', fgz);
        FileExists = true;
        f = fgz;
    else
        FileExists = false;
    end
end
end

%==========================================================================
function G = updateG(G,v,labels)
I = G.indices;
G.Aind_gm_ctx_lh         = ismember(find(v==3), I);
G.Aind_gm_ctx_rh         = ismember(find(v==42), I);
G.Aind_gm_ctx_unassigned = ismember(find(v==50), I);
G.indices_gm_ctx_lh = G.indices(G.Aind_gm_ctx_lh);
G.indices_gm_ctx_rh = G.indices(G.Aind_gm_ctx_rh);
if ismember(50, labels)
    G.indices_gm_ctx_unassigned = G.indices(G.Aind_gm_ctx_unassigned);
else
    G.indices_gm_ctx_unassigned = [];
end
end

%==========================================================================
function v_rib = getrib(G, f_gmctx, DirSubjects)
h_gmctx = spm_vol(f_gmctx);
f_rib = G.f.source1;
if not(isempty(DirSubjects))
    if isfield(G.f, 'hcpsave_root')
        f_hcp = G.f.hcpsave_root;
    else
        f_hcp = G.f.hcp_root;
    end
    f_rib = strrep(f_rib, f_hcp, DirSubjects);
end
f_tmp = strrep(f_rib, '.nii', '_tmp.nii');
hb_nii_reslice(f_rib, f_gmctx, 0, f_tmp);
h_tmp = spm_vol(f_tmp);
assert(isequal(h_tmp.dim,h_gmctx.dim));
assert(all(abs(h_tmp.mat-h_gmctx.mat)<1e-6,'all'));
v_rib = spm_read_vols(h_tmp);
delete(f_tmp);
end

%==========================================================================
function v = getlhrhvol(v_gmctx, v_rib, Ngmctx)
v = zeros(size(v_gmctx));
se = strel('sphere',1);

% lh
vlh = v;
vlh(v_rib==3) = 1;
vlh = imdilate(vlh, se);
vlh(and(v_gmctx, logical(vlh))) = 3;

% rh
vrh = v;
vrh(v_rib==42) = 1;
vrh = imdilate(vrh, se);
vrh(and(v_gmctx, logical(vrh))) = 42;

% handle common labels
cmn = and(logical(vlh), logical(vrh));
assert((nnz(cmn)/Ngmctx)<0.02, 'fishy: too many common labels'); % 2% is reasonable
vrh(cmn) = 0; % label common ones as lh
v(vlh==3) = 3;
v(vrh==42) = 42;
end

%==========================================================================
function UniqueLabels = writelhrhvol(v, v_gmctx, Ngmctx, h_gmctx, f_o)
if nnz(v)==Ngmctx
    % all good
    UniqueLabels = [0 3 42]; % [non-mask lh rh]
else
    Ndiff = Ngmctx - nnz(v);
    assert(Ndiff>0, 'fishy: nnz(v) should be < Ngmctx');
    assert(Ndiff<10, 'fishy: too many voxels missed out');
    vdiff = logical(v_gmctx)-logical(v);
    I = find(vdiff);
    if 1
        v(I) = 50; % fishy voxels labelled as 50
        UniqueLabels = [0 3 42 50]; % [non-mask lh unknown rh]
    else
        %-this approach is not stable--------------------------------------
        for k=1:Ndiff
            sesz = 1;
            while sesz>0
                v0 = zeros(size(v));
                v0(I(k)) = 1;
                v0 = imdilate(v0, strel('sphere',sesz));
                nlh = nnz(and(logical(v0), v==3));
                nrh = nnz(and(logical(v0), v==42));
                if nlh==nrh
                    sesz = sesz+1;
                else
                    if nlh>nrh % more lh neighbours
                        v(I(k)) = 3;
                    else  % more rh neighbours
                        v(I(k)) = 42;
                    end
                    sesz = 0;
                end
            end
        end
        UniqueLabels = [0 3 42 50];
    end
    assert(nnz(v)==Ngmctx);
end

assert(isequal(sort(unique(v(:))),UniqueLabels(:)));

h = struct;
h.fname = f_o;
h.dim = h_gmctx.dim;
h.mat = h_gmctx.mat;
h.dt = [16 0]; % [2 0] does not prob´perly save values as integers
spm_write_vol(h,v);
gzip(f_o);
delete(f_o);

end

%==========================================================================
function [h_gmctx, v_gmctx] = getgmctx(G, f_gmctx, DirSubjects)
Ngmctx = length(G.indices_gm_ctx);
if exist(f_gmctx, 'file')
    GunzipDlt = false;
else
    if exist([f_gmctx,'.gz'], 'file')
        ThrowErr = false;
        gunzip(f_gmctx);
        GunzipDlt = true;
    else
        GunzipDlt = false;
        if isempty(DirSubjects)
            ThrowErr = true;
        else
            if isfield(G.f, 'hcpsave_root')
                f_hcp = G.f.hcpsave_root;
            else
                f_hcp = G.f.hcp_root;
            end
            d =  strrep(f_gmctx, f_hcp, DirSubjects);
            if exist(d, 'file')
                f_gmctx = d;
            else
                ThrowErr = true;
            end
        end
    end
    if ThrowErr
        error('Missing file: %s',f_gmctx);
    end
end

h_gmctx = spm_vol(f_gmctx);
v_gmctx = spm_read_vols(h_gmctx);
assert(nnz(v_gmctx)==Ngmctx);

if GunzipDlt
    assert(exist([f_gmctx, '.gz'],'file'));
    delete(f_gmctx);
end
end