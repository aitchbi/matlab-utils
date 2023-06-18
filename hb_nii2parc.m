function [p, E, l, P] = hb_nii2parc(f_i, f_p, varargin)
% HB_NII2PARC parcellates a given nifti based on a given parcellation.
%
% 
% Inputs:
%   f_i: a nifti volume; current code only works for 3D volume, not 4D.
% 
%   f_p: a segmentation/parcellation nifti file, e.g. a mask or a
%   cortical parcellation.
%
%   Name-Value Pair Arguments:
%   RegisterToParc: logical [default: false]. By default, the parcellation
%   file f_p will be registered (affine transformed) to f_i so that they
%   match voxel-to-voxel, in case this already is not the case. However,
%   the operation can be reversed, that is, to register the parcellation
%   file to the f_i.  
%
%   DeleteIntermediateFiles: logical [default: false]. 
%
%   InterpolationOrder: 0 or 1. If not set, y default, 0 used if f_p needs
%   reslicing whereas 1 is used if f_i needs reslicing. This is because f_p
%   is a segmentation, and we don't want to mess up the label, thus, we use
%   nearest neighbor interpolation, whereas for f_i, typically that's not a
%   segmentation, thus, values are continious, and therefore, best to use
%   linear interpolation. But if f_i is a segmentation (label, binary
%   file), set this input to 0.
%
%   MaskFile: filepath to a nifti mask file, i.e., binary. f_i will be
%   masked, before parcellating it with f_p. The mask can be either in
%   register with f_i, f_p, or even not in register with either; the
%   function reslices the mask appropriately if needed.
%
% Outputs:
%   p: fraction of total energy of f_i in each parcel of f_p. If no mask
%   file specified: first, the total energy (norm^2) of f_i across all
%   parcels in f_p is computed, denoted E (see below, second output), then,
%   the enrgy associated to each parcel is computed and returned, such that
%   E = p(1) + p(2) + ... + p(Np), where Np denotes the number of parcels.
%   If mask file specified, the procedure is the same as above, only that
%   initially f_p is masked by the mask file; the equality condition still
%   holds.
%
%   E: total energy (norm^2) of f_i across all parcels in f_p, after
%   registereing f_i and f_p, and potentailly applied an input mask.
%   Specifically, TE is the sum of second power of all voxels in f_i that
%   overlap with a parcel in f_p and are also within the mask if given.
%
%   l: the labels in f_p associated to each element of P; the labels are
%   sorted in ascending order, i.e, L(i) < L(i+1) for all i \in [1, Np]
%   where Np denotes the number of unique non-zero labels in f_p.
%
%   P: an MxF cell array, wherein array (i,j) stores voxel values of j-th
%   frame of f_i associated to label l(j) in f_p. The array is sorted based
%   on sorted labels in f_p; i.e., P{1} and P{end} contain voxels in f_i
%   associated to the smallest and largest non-zero label in f_p,
%   respectively.
% 
% Dependencies:
% .SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_nii_reslice.m 
% .https://github.com/aitchbi/matlab-utils/blob/main/spm_modified.m
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_nii_handlegzip.m
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_nii_load.m 
% 
% Hamid Behjat

%-process optional inputs--------------------------------------------------
funcLogi = @(x) assert(islogical(x));
funcIntp = @(x) assert(ismember(x,[0 1]) || isempty(x));
funcPath = @(x) assert(ischar(x) || isempty(x));
p = inputParser;
addParameter(p,'RegisterToParc', false, funcLogi);
addParameter(p,'DeleteIntermediateFiles', true, funcLogi);
addParameter(p,'InterpolationOrder', [], funcIntp);
addParameter(p,'MaskFile', [], funcPath);
parse(p,varargin{:});
opts = p.Results;
%--------------------------------------------------------------------------

[f_i, f_p, FilesToCleanUp] = verifyinputfile(f_i, f_p);

%assert(length(spm_vol(f_i))==1, 'extend for 4D');

%-Verify inputs & register if needed.
[f_i, f_p, FilesToCleanUp] = registerfiles(opts, f_i, f_p, FilesToCleanUp);

%-Mask files if mask given.
[f_i, f_p, FilesToCleanUp, I_msk] = maskfiles(opts, f_i, f_p, FilesToCleanUp);

%-Load files.
v_i = hb_nii_load(f_i, 'IndicesToLoad', I_msk);
v_p = hb_nii_load(f_p, 'IndicesToLoad', I_msk);

%-Initialize outputs, etc. 
if length(size(v_i))==4
    Nv = size(v_i,4);
else
    Nv = 1;
end
l = sort(unique(v_p)); % labels
assert(l(1)==0, 'fishy');
l = l(2:end);
Np = length(l);
P = cell(Np,Nv);
E = zeros(1,Nv);
p = zeros(Np,Nv);

for iV=1:Nv

    %-Extract a volume.
    if Nv==1
        viv = v_i;
    else
        viv = v_i(:,:,:,iV);
        showprgs(iV,Nv,'Parcellating volumes..');
    end

    %-Compute total energy.
    vivp = viv(v_p~=0);
    E(iV) = norm(vivp)^2;
    
    %-Parcellate & compute energies.
    for iP=1:Np
        P{iP,iV} = viv(v_p==l(iP));  
        p(iP,iV) = norm(P{iP,iV})^2;
    end
    assert(abs(sum(p(:,iV))-E(iV))<1e-6*E(iV));
    assert(sum(cellfun(@length, P(:,iV)))==length(vivp));
end

%-Cleanup.
cleanup(opts, FilesToCleanUp);
end

%==========================================================================
function [f_i, f_p, FilesToCleanUp] = verifyinputfile(f_i, f_p)
[f_i, f_dlt1, sts1] = hb_nii_handlegzip(f_i);
[f_p, f_dlt2, sts2] = hb_nii_handlegzip(f_p);
assert(sts1==1, 'Problem with input file f_i [sts: %d]', sts1);
assert(sts2==1, 'Problem with input file f_p [sts: %d]', sts2);
FilesToCleanUp = [];
FilesToCleanUp = appendftc(FilesToCleanUp, f_dlt1);
FilesToCleanUp = appendftc(FilesToCleanUp, f_dlt2);
end

%==========================================================================
function sts = chkmatch(f1,f2)
tol = 1e-6;
h1 = spm_vol(f1);
h2 = spm_vol(f2);
if length(h1)>1
    h1 = h1(1);
end
if length(h2)>1
    h2 = h2(1);
end
d1 = isequal(h1.dim, h2.dim);
d2 = all(all(abs(h1.mat-h2.mat)<tol));
if d1 && d2
    sts = true;
else
    sts = false;
end
end

%==========================================================================
function [f_i, f_p, FilesToCleanUp] = registerfiles(opts, f_i, f_p, FilesToCleanUp)

sts = chkmatch(f_i, f_p);

if sts
    % f_i & f_p are matched voxel-to-voxel
else
    % reslice f_i/f_p so they match voxel-to-voxel
    
    if opts.RegisterToParc % transform f_i
        
        if isempty(opts.InterpolationOrder)
            ord = 1; % linear interp
        else
            ord = opts.InterpolationOrder;
        end
        f_reg = strrep(f_i, '.nii', '_tmp.nii');
        hb_nii_reslice(f_i, f_p, ord, f_reg);
        assert(chkmatch(f_reg, f_p));
        f_i = f_reg;
    
    else % transform f_p
        if isempty(opts.InterpolationOrder)
            ord = 0; % nearest neighb interp, to not mess up labels
        else
            ord = opts.InterpolationOrder;
        end
        f_reg = strrep(f_p, '.nii', '_tmp.nii');
        hb_nii_reslice(f_p, f_i, ord, f_reg);
        assert(chkmatch(f_reg, f_i));
        f_p = f_reg;
    end
    FilesToCleanUp = appendftc(FilesToCleanUp, f_reg);
end
end

%==========================================================================
function F = appendftc(F, f)
if isempty(f)
    return;
end
if isempty(F)
    F = {f};
else
    F = [
        F
        f
        ];
end
end

%==========================================================================
function [f_im, f_pm, FTC, I_msk] = maskfiles(opts, f_i, f_p, FTC)
if isempty(opts.MaskFile)
    f_im  = f_i;
    f_pm  = f_p;
    I_msk = [];
else
    [f_msk, f_dlt, sts] = hb_nii_handlegzip(opts.MaskFile);
    assert(sts==1, 'Problem with mask file [sts: %d]', sts);
    FTC = appendftc(FTC, f_dlt);
    sts = chkmatch(f_msk, f_p);
    if sts
        % mask already in register with f_p & f_i
        f_mskr = f_msk;
    else
        % reslice mask to match f_p (or f_i, as both should be in match)
        f_mskr = strrep(f_msk, '.nii', '_tmp.nii');
        hb_nii_reslice(f_msk, f_p, 0, f_mskr);
        assert(chkmatch(f_mskr, f_i));
        assert(chkmatch(f_mskr, f_p));
        FTC = appendftc(FTC, f_mskr);
    end
    f_im = maskfiles_pvt(f_i, f_mskr);
    f_pm = maskfiles_pvt(f_p, f_mskr);
    FTC = appendftc(FTC, f_im);
    FTC = appendftc(FTC, f_pm);
    I_msk = find(spm_read_vols(spm_vol(f_mskr)));
end
end

%==========================================================================
function f_im = maskfiles_pvt(f_i, f_m)
assert(logical(chkmatch(f_i, f_m)));
f_im = strrep(f_i,'.nii', '_masked.nii');
v_m = spm_read_vols(spm_vol(f_m));
I_m = find(v_m);
[v_i, h_i] = hb_nii_load(f_i,'IndicesToLoad',I_m);    
Nv = length(h_i);
if Nv==1
    h_ref = h_i;
    v0 = zeros(size(v_i));
else
    h_ref = h_i(1);
    v0 = zeros(size(v_i(:,:,:,1)));
end
for iV=1:Nv
    if Nv==1
        viv = v_i;
    else
        viv = v_i(:,:,:,iV);
        showprgs(iV,Nv,'Masking volumes..');
    end
    v_im = v0;
    v_im(I_m) = viv(I_m);
    if iV==1
        h_im = struct;
        h_im.fname = f_im;
        h_im.dim = h_ref.dim;
        h_im.mat = h_ref.mat;
        h_im.dt  = h_ref.dt;
        h_im = spm_create_vol(h_im);
    end
    h_im.n(1) = iV;
    spm_write_vol(h_im, v_im);
end
end

%==========================================================================
function cleanup(opts,F)
if opts.DeleteIntermediateFiles
    if not(isempty(F))
        for k=1:length(F)
        delete(F{k});
        end
    end
end
end

%==========================================================================
function showprgs(n,N,tag)
l = numel(num2str(N));
if n==1
    fprintf('\n..%s ',tag);
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end
