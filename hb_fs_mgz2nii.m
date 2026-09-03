function hb_fs_mgz2nii(f_mgz, f_nii, dir_subjs, dir_freesurfer, dir_hbfssh)
%
%
%
% f_orig eg: subjdir/mri/ribbon.mgz
%
% f_gii eg: subjdir/mri/ribbon.nii
%
% dir_subjs: directory where FS subject folders are in; eg:
%   on gonzo: .../bfdata_working_copy/fs/x
% on aitchbi: .../HCP_database
%
% dir_freesurfer: where FS is installed; eg:
%   on gonzo: '/opt/fs72' 
% on aitchbi: '/Applications/freesurfer/7.1.1'
%
% dir_hbfssh: directory where hb_fs_mris_convert.sh is located
%
% h behjat

assert(endsWith(f_mgz, '.mgz'), 'source file not mgz');

assert(endsWith(f_nii, '.nii'), 'destination file not nii');

assert(exist(f_mgz, 'file'));

cmd = sprintf('%s -d %s -r %s -i %s -o %s', ...
    fullfile(dir_hbfssh, 'hb_fs_mri_convert.sh'), ...
    dir_subjs, ...
    dir_freesurfer, ...
    f_mgz, ...
    f_nii ...
    );
hb_runcmd(cmd, 'mgz to nii conversion failed.');
end