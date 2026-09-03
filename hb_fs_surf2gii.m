function hb_fs_surf2gii(f_orig, f_gii, dir_subjs, dir_freesurfer, dir_hbfssh)
%
%
% h behjat 
%
%
% f_orig eg:
% read-dir/surf/lh.pial
% read-dir/surf/rh.white

% f_gii eg:
% write-dir/surf/lh.pial.surf.gii
% write-dir/surf/rh.white.surf.gii
%
% dir_subjs: directory where FS subject folders are in; eg:
%   on gonzo: .../bfdata_working_copy/fs/x for BF
% on aitchbi: .../HCP_database
%
% dir_freesurfer: where FS is installed; eg:
%   on gonzo: '/opt/fs72' 
% on aitchbi: '/Applications/freesurfer/7.1.1'
%
% dir_hbfssh: directory where hb_fs_mris_convert.sh is located
%
%
% h behjat

assert(exist(f_orig, 'file'));

assert(endsWith(f_gii, '.surf.gii'));

cmd = sprintf('%s -d %s -r %s -i %s -o %s',...
    fullfile(dir_hbfssh, 'hb_fs_mris_convert.sh'),...
    dir_subjs,...
    dir_freesurfer,...
    f_orig,... % xx.pial/white
    f_gii...   % xx.pial/white.surf.gii
    );

hb_runcmd(cmd,'Error converting surface file.');
end