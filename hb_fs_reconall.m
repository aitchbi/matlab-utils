function hb_fs_reconall(f_t1, subj_name, dir_subjs, dir_freesurfer, dir_hbfssh)
%
% inputs: 
%   f_t1: full path to t1w image (*.nii or *.nii.gz)
%
%   subj_name: subject name. freesurfer outputs will bve saved in:
%       fullfile(dir_subjs, subj_name)
%
%   dir_subjs: directory where FS subject folders are in; eg:
%       on gonzo: .../bfdata_working_copy/fs/x for BF
%       on aitchbi: .../HCP_database
%
%   dir_freesurfer: where FS is installed; eg:
%       on gonzo: '/opt/fs72' 
%       on aitchbi: '/Applications/freesurfer/7.1.1'
%
% dir_hbfssh: directory where hb_fs_mris_convert.sh is located
%
%
% h behjat

assert(exist(f_t1, 'file'));
assert(ischar(subj_name));

fprintf('\n.running freesurfer recon-all on: %s [started: %s]\n', ...
    f_t1, ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));

cmd = sprintf('%s -d %s -r %s -i %s -s %s',...
    fullfile(dir_hbfssh, 'hb_fs_reconall.sh'),...
    dir_subjs,...
    dir_freesurfer,...
    f_t1, ...
    subj_name...
    );

t = hb_runcmd(cmd, 'error running freesurfer recon-all');

d = 'recon-all successful | total time:';
if t < 180       % < 3 minutes
    fprintf('\n.%s %d seconds\n', d, round(t));
elseif t < 10800 % < 3 hours
    fprintf('\n.%s %d minutes\n', d, round(t/60));
else             % > 3 hours
    fprintf('\n.%s %0.1f hours\n', d, t/3600);
end
end