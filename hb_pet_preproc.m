function [f_pet3Dt1, sts] = hb_pet_preproc(f_pet4D, f_t1, varargin)
% hb_pet_preproc
%
% inputs
%   f_pet4D: absolute path to 4D tau-PET nifti (moved)
%   f_t1   : absolute path path to T1w nifti (not moved)
%
% output
%   f_pet3Dt1: 4D PET motion-corrected, then summed to get a 3D,
%   then registered & resliced to match input T1w space/grid
%   sts: 1 (file generated) or 2 (file exists).
%
% dependencies
%   SPM12
%   github/aitchbi/matlab-utils/hb_nii_reslice_v2.m
%   github/aitchbi/matlab-utils/misc/spm_modified/spm_run_coreg_hb.m
% 
% h behjat

p = inputParser;
addParameter(p, 'DeleteAuxiliaryFiles', true, @(x) islogical(x));
addParameter(p, 'SaveSanityCheckFigures', false, @(x) islogical(x));
addParameter(p, 'DebugMode', false, @(x) islogical(x));
addParameter(p, 'JustGetOutputFilename', false, @(x) islogical(x));
parse(p,varargin{:});
opts = p.Results;

if endsWith(f_pet4D, '.gz')
    gunzip(f_pet4D);
    f_pet4D = strrep(f_pet4D, '.gz', '');
    ToDelete = f_pet4D;
else
    ToDelete = [];
end

% step 0: names------------------------------------------------------------
% names:
[d_pet, n_pet]          = fileparts(f_pet4D);
[d_t1, n_t1, e_t1]      = fileparts(f_t1);

n_t1_recentered         = [n_t1, '-recentered'];
n_pet_mc_native_v0      = ['motioncorr_', n_pet];
n_pet_mc_native         = [n_pet, '.motioncorr'];
n_pet_mc_native_sum     = [n_pet_mc_native, '.sum'];

tag_rt1                 = sprintf('.realigned2-%s-orig', n_t1);
tag_rt1tmp              = sprintf('.realigned2-%s-recentered', n_t1);
tag_rt1tmp_prefix       = sprintf('realigned2-%s-recentered_', n_t1);

f_petmat                = strrep(f_pet4D, '.nii', '.mat'); % ? (an auxillary output from step 1)

f_t1_recentered         = fullfile(d_t1, [n_t1_recentered, e_t1]);

f_pet_rp_txt_v0         = fullfile(d_pet, ['rp_', n_pet, '.txt']); % realignment params
f_pet_rp_txt            = fullfile(d_pet, [n_pet, '.realignmentParams.txt']); 

f_pet_mc_native_v0      = fullfile(d_pet, [n_pet_mc_native_v0, '.nii']); % realigned (motion-corrected) 4D PET
f_pet_mc_native         = fullfile(d_pet, [n_pet_mc_native, '.nii']); 

f_pet_mc_sum_native     = fullfile(d_pet, [n_pet_mc_native_sum, '.nii']);

f_pet_mc_mean_native_v0 = fullfile(d_pet, ['mean', n_pet, '.nii']); % output name format from SPM
f_pet_mc_mean_native    = fullfile(d_pet, [n_pet, '.mean.nii']); % mean PET output from realignment

f_pet_mc_sum_inT1tmp_v0 = fullfile(d_pet, [tag_rt1tmp_prefix, n_pet_mc_native_sum, '.nii']); % resliced to match tmp T1 (T1 recentered to center of mass)
f_pet_mc_sum_inT1tmp    = fullfile(d_pet, [n_pet_mc_native_sum, tag_rt1tmp, '.nii']);

f_pet_mc_sum_inT1       = fullfile(d_pet, [n_pet_mc_native_sum, tag_rt1, '.nii']); % resliced to match original input T1 

d_snap                  = fullfile(d_pet, 'snapshots');
f_png1                  = fullfile(d_snap, 'sanity-check-1_pet-and-t1-centered-to-cent-of-mass.png');
f_png2                  = fullfile(d_snap, 'sanity-check-2_pet-and-t1-registred-resliced-to-match-recentered-t1-space.png');
f_png3                  = fullfile(d_snap, 'sanity-check-3_pet-and-t1-registred-resliced-to-match-input-t1-space.png');
f_pet3Dt1               = [f_pet_mc_sum_inT1, '.gz']; % output

chk1 = exist(f_pet3Dt1, 'file');
chk2 = exist(f_pet_mc_sum_inT1, 'file');

if chk1
    sts = 2;
elseif chk2
    gzip(f_pet_mc_sum_inT1);
    delete(f_pet_mc_sum_inT1);
    sts = 2;
end

if opts.JustGetOutputFilename
    if ~exist('sts', 'var')
        sts = [];
    end
    return;
end

if chk1 || chk2
    fprintf('\n.preprocessed pet file exists: %s\n', f_pet3Dt1);
    return;
end

[~,~] = mkdir(d_snap);

copyfile(f_t1, f_t1_recentered);

% init SPM ----------------------------------------------------------------
spm('Defaults','PET');
spm_jobman('initcfg');

% expand 4D frames --------------------------------------------------------
h_pet4D = spm_vol(f_pet4D);
nFrames = numel(h_pet4D);
if nFrames < 2
    error('PET must be a 4D nifti; found %d volume(s).', nFrames);
end

petFrames = cell(nFrames,1);
for i=1:nFrames
    petFrames{i} = sprintf('%s,%d', f_pet4D, i);
end

% step 1: realign (motion correction) -------------------------------------
matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.estwrite.data = {petFrames};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 1; % align to mean
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which   = [2 1]; % reslice all + mean
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp  = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'motioncorr_';
spm_jobman('run', matlabbatch);
% input: <name>.nii
% outputs: 
%   <name>.mat      : ?
%   rp_<name>.txt   : f_rp_txt
%   r<name>.nii     : f_pet_mc_native
%   mean<name>.nii  : f_pet_mc_native_mean
movefile(f_pet_rp_txt_v0, f_pet_rp_txt);
movefile(f_pet_mc_native_v0, f_pet_mc_native);
movefile(f_pet_mc_mean_native_v0, f_pet_mc_mean_native);

% step 2: sum realigned frames to a 3D volume (native PET space) ----------
rFrames = cell(nFrames, 1);
for i=1:nFrames
    rFrames{i} = fullfile(d_pet, [n_pet_mc_native, '.nii,' num2str(i)]);
end
h_pet_mc = spm_vol(char(rFrames));
v_pet_mc = spm_read_vols(h_pet_mc);
v_sum = sum(v_pet_mc, 4);
h_pet_mc_sum = h_pet_mc(1);
h_pet_mc_sum.fname = f_pet_mc_sum_native;
h_pet_mc_sum.dt    = [spm_type('float32') 0];
spm_write_vol(h_pet_mc_sum, v_sum);

% step 3A: center images---------------------------------------------------
% center images so there 0 0 0 coordinate is close
% otherwise registration via SPM fails if images poorly overlap
% the PET and T1 should be close in field of view in world coordinates
recenter_to_center_of_mass(f_t1_recentered);
recenter_to_center_of_mass(f_pet_mc_mean_native);
recenter_to_center_of_mass(f_pet_mc_sum_native);
if opts.DebugMode || opts.SaveSanityCheckFigures
    spm_check_registration( ...
        f_t1_recentered, ...
        f_pet_mc_mean_native, ...
        f_pet_mc_sum_native); 
    % should "roughly" overlap
    if opts.SaveSanityCheckFigures
        exportgraphics(gcf, f_png1, 'Resolution', 300)
        close(gcf);
    end
end

% step 3B: coregister: estimate only --------------------------------------
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.estimate.ref    = {f_t1_recentered};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {f_pet_mc_mean_native};
matlabbatch{1}.spm.spatial.coreg.estimate.other  = {f_pet_mc_sum_native};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [
    0.02  * ones(3,1) % x,y,z translations [mm]
    0.001 * ones(3,1) % rotations (pitch, roll, yaw) [radians]
    0.01  * ones(3,1) % zooms [unitless]
    0.001 * ones(3,1) % shears [unitless]
    ]';
spm_jobman('run', matlabbatch);

% step 3C: reslice summed PET onto mass-centered T1 grid ------------------
P = char(f_t1_recentered, f_pet_mc_sum_native);
flags = struct();
flags.mask   = 0;
flags.mean   = 0;
flags.interp = 4;
flags.which  = 1; % reslice second image only
flags.wrap   = [0 0 0];
flags.prefix = tag_rt1tmp_prefix;
spm_reslice(P, flags);
movefile(f_pet_mc_sum_inT1tmp_v0, f_pet_mc_sum_inT1tmp);

if opts.DebugMode || opts.SaveSanityCheckFigures
    spm_check_registration( ...
        f_t1_recentered, ...
        f_pet_mc_sum_inT1tmp); 
    % should "closely" overlap [in input re-centered T1w space/grid]
    if opts.SaveSanityCheckFigures
        exportgraphics(gcf, f_png2, 'Resolution', 300)
        close(gcf);
    end
end

% step 3D: reslice output from 3C onto "orig input T1" grid ---------------
%------3D-1: get reg matrix for: f_t1_tmp → f_t1 
f_i = f_t1_recentered;
f_r = f_t1;
f_o = strrep(f_t1_recentered, '.nii', sprintf('_tmp%d.nii', randi([1e5, 999999])));
[~, M] = hb_nii_reslice_v2(f_i, f_r, f_o, 'RegisterThenReslice', true);
if opts.DebugMode
    spm_check_registration(f_t1, f_o); % should "closely" overlap
end
delete(f_o);

%------3D-2: use reg matrix to: f_pet_mc_sum_inT1tmp → f_t1
f_i = f_pet_mc_sum_inT1tmp;
f_r = f_t1;
f_o = f_pet_mc_sum_inT1;
hb_nii_reslice_v2(f_i, f_r, f_o, ...
    'RegisterThenReslice', true, ...
    'RegisterationMatrix', M);
if opts.DebugMode || opts.SaveSanityCheckFigures
    spm_check_registration( ...
        f_t1, ...
        f_pet_mc_sum_inT1); 
    % should "closely" overlap [in input T1w space/grid]
    if opts.SaveSanityCheckFigures
        exportgraphics(gcf, f_png3, 'Resolution', 300)
        close(gcf);
    end
end

% clean up-----------------------------------------------------------------
F = {
    f_petmat             % ? an auxiliary output
    f_pet_rp_txt         % realignment parameters file (if produced)
    f_pet_mc_native      % realigned 4D PET in native space
    f_pet_mc_mean_native % mean PET in native space from realignment
    f_pet_mc_sum_inT1tmp % sum PET resliced to recentered t1w space
    f_t1_recentered      % input T1 recentered to its center of mass (tmp file)
    };
if opts.DeleteAuxiliaryFiles
    for k=1:length(F)
        delete(F{k});
    end
else
    for k=1:length(F)
        fk = F{k};
        if endsWith(fk, '.nii')
            gzip(fk);
            delete(fk);
        end
    end
end

% to always keep:
% f_pet_mc_sum_inT1   % motion-corrected summed PET registered & resliced to input T1 space
% f_pet_mc_sum_native % motion-corrected summed PET in native PET space

gzip(f_pet_mc_sum_native); 
delete(f_pet_mc_sum_native);

gzip(f_pet_mc_sum_inT1);
delete(f_pet_mc_sum_inT1);

if ~isempty(ToDelete)
    delete(ToDelete);
end

sts = 1;
end

%==========================================================================
function recenter_to_center_of_mass(fname)
    h = spm_vol(fname);
    v = spm_read_vols(h);

    % compute center of mass (non-zero voxels)
    [x,y,z] = ind2sub(size(v), find(v > 0));
    com = [mean(x), mean(y), mean(z), 1]';

    % voxel -> world
    com_world = h.mat * com;

    % shift origin to center of mass
    hmat = h.mat;
    hmat(1:3,4) = hmat(1:3,4) - com_world(1:3);

    spm_get_space(fname, hmat);
end
