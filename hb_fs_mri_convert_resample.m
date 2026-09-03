function hb_fs_mri_convert_resample(f_i, f_o, voxsize, iorder, dir_subjs, dir_freesurfer, dir_hbfssh)
%
%
%
% f_i: eg subjdir/mri/T1.nii
%
% f_o: eg subjdir/mri/T1_2mm.nii
%
% voxsize: eg 2
% 
% iorder: interpolation order; 1 or 2
%
% dir_subjs: directory where FS subject folders are in; eg:
%   on gonzo: .../bfdata_working_copy/fs/x
% on aitchbi: .../HCP_database
%
% dir_freesurfer: where FS is installed; eg:
%   on gonzo: '/opt/fs72' 
% on aitchbi: '/Applications/freesurfer/7.1.1'
%
% dir_hbfssh: directory where hb_fs_mri_convert_resample.sh is located
%
% h behjat

assert(exist(f_i, 'file'));
assert(isnumeric(voxsize));
assert(length(voxsize)==1, 'only isotropic voxels'); % NOTE 1
assert(ismember(iorder, [0 1]));

cmd = sprintf(...
    '%s -d %s -r %s -i %s -o %s -v %d -t %d',...
    fullfile(dir_hbfssh, 'hb_fs_mri_convert_resample.sh'),...
    dir_subjs,...
    dir_freesurfer,...
    f_i,...
    f_o,...
    voxsize, ...
    iorder); % isotropic vox size

hb_runcmd(cmd, 'resampling via mri_convert failed.');
end

% NOTE 1
% hb_fs_mri_convert_resample.sh only takes a scalar as input for voxel
% resolution; need to extend to allow to take as input a 1x3 vector, and
% then use a switch in the code for either scenario.
