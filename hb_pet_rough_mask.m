function hb_pet_rough_mask(f_pet_in, f_pet_out, varargin)
% rough_pet_mask_nifti
%
% inputs:
%   f_pet_in  : absolute path to input PET nifti (.nii or .nii.gz)
%   f_pet_out : absolute path to output masked PET nifti
%
% optional Name-Value:
%   'WriteMask' : true/false (default: false)
%   'MaskFilename': path to output mask 
%                   (default: strrep(output_pet, '.nii', '_mask.nii'))
%
% dependencies:
%   SPM12
%
% h behjat 

p = inputParser;
addParameter(p, 'WriteMask', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'FilenameMask', [], @(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opts = p.Results;

input_pet_nii = ensure_nii_unzipped(f_pet_in);

h_in = spm_vol(input_pet_nii);
v_in = spm_read_vols(h_in);
v_in = double(v_in);

%-compute threshold.
vals = v_in(v_in > 0);
if isempty(vals)
    error('PET image appears empty');
end
thr = 0.1 * prctile(vals, 95);
mask = v_in > thr;

%-keep largest component.
CC = bwconncomp(mask, 6);
if CC.NumObjects > 0
    nvox = cellfun(@numel, CC.PixelIdxList);
    [~, idx] = max(nvox);
    mask2 = false(size(mask));
    mask2(CC.PixelIdxList{idx}) = true;
    mask = mask2;
end

%-fill holes slice-wise.
for z = 1:size(mask,3)
    mask(:,:,z) = imfill(mask(:,:,z), 'holes');
end
mask = imdilate(mask, strel('sphere',1)); % optional dilation
v_in(~mask) = 0; % apply mask

%-write output.
h_out = h_in;
h_out.fname = strip_gz(f_pet_out); % write as .nii first
h_out.dt = [spm_type('float32') 0];
spm_write_vol(h_out, v_in);

%-gzip.
if endsWith(f_pet_out, '.gz')
    gzip(h_out.fname);
    delete(h_out.fname);
end

%-write mask.
if opts.WriteMask
    h_mask = struct;
    h_mask.dim = h_in.dim;
    h_mask.mat = h_in.mat;
    h_mask.dt  = [2 0];
    [pdir, pname, ~] = fileparts(strip_gz(f_pet_out));
    if isempty(opts.FilenameMask)
        f_mask = fullfile(pdir, [pname '_mask.nii']);
    else
        f_mask = opts.FilenameMask;
    end
    h_mask.fname = strip_gz(f_mask);
    spm_write_vol(h_mask, double(mask));
    if endsWith(f_mask, '.gz')
        gzip(h_mask.fname);
        delete(h_mask.fname);
    end
end
end

%==========================================================================
function fnii = ensure_nii_unzipped(f)
    f = char(f);
    if endsWith(f, '.nii.gz')
        fnii = f(1:end-3);
        if exist(fnii,'file') ~= 2
            gunzip(f);
        end
    else
        fnii = f;
    end
end

%==========================================================================
function f = strip_gz(f)
    if endsWith(f, '.nii.gz')
        f = f(1:end-3);
    end
end