function [f_surfs, Y] = hb_vol2surf(f_vol,hemi,fsopts,varargin)
% HB_VOL2SURF projects a given nifti volume to surface. 
%
% surface projection can be done to: 
%  - subject's own surface (uisng FreeSurfer)
%  - fsaverage surface     (using FreeSurfer)
%  - fs_LR surface         (using Connectome Workbench)
%
% note: when mapping to fs_LR, first, mapping is done to fsaverage, and
% then from fsaverage to fs_LR.
%
% note: projection is done to white surface; though i guess projection to
% e.g. the pial etc. is not different, just a difference in visualisation.
%
% HB

d = inputParser;
addParameter(d, 'SurfaceSpace', 'subject');
addParameter(d, 'FramesToMap', []);
addParameter(d, 'SaveSurface', true);
addParameter(d, 'ReturnAsVector', false);
addParameter(d, 'OutputSurfaceName', []);
addParameter(d, 'Silent', true);
parse(d, varargin{:});
opts = d.Results;
    
[f_vol, hemis, f_surfs, F_cleanup] = chkinputs(f_vol, hemi, fsopts, opts);

for iH=1:length(hemis)

    hemi = hemis{iH};

    f_surf = f_surfs.(hemi);

    switch opts.SurfaceSpace
        case {'subject', 'fsaverage'}
            f_surfstp1 = f_surf;
        case 'fsLR'
            f_surfstp1 = strrep(f_surf, 'fsLR.gii', 'fsaverage.gii');
    end

    cmd = sprintf(...
        '%s -s %s -i %s -h %s -b %s -r %s -o %s',...
        fsopts.sh_vol2surf,...
        fsopts.ID,...              %s
        f_vol,...                  %i
        hemi,...                   %h
        fsopts.dir_subjs,...       %b
        fsopts.dir_freesurfer, ... %r
        f_surfstp1);               %o

    switch opts.SurfaceSpace
        case {'fsaverage', 'fsLR'}
            cmd = sprintf('%s -d %s', cmd, 'fsaverage_hb');
    end
    % e.g. mri_vol2surf --srcsubject 100307 --trgsubject 100307 --hemi lh --mov 'test.nii' --regheader 100307 --o 'lh_test.gii' --projfrac 0.5
    
    runcmd(cmd, 'FreeSurfer mri_vol2surf', f_surfstp1);
    
    switch opts.SurfaceSpace
        case {'subject', 'fsaverage'}
            continue;
            
        case 'fs_LR'
            % transform map in fsaverage to fs_LR
            % f_surfstp1 -> f_surf
            cmd = sprintf('neuromaps.transforms.fsaverage_to_fslr(%s, target_density=''32k'', hemi=None, method=''linear''))', f_surfstp1);
            runcmd(cmd, 'Neuromaps transforms.fsaverage_to_fslr', []);
            %pyrun(cmd);
    end
end

docleanup(F_cleanup)

if opts.ReturnAsVector
    error('code..');
else
    Y = [];
end
end

%==========================================================================
function [f_vol,hemis,f_surfs,F_cleanup] = chkinputs(f_vol,hemi,fsopts, opts)

F_cleanup = [];

assert(endsWith(f_vol, {'.nii', '.nii.gz'}));

if endsWith(f_vol, '.nii.gz')
    gunzip(f_vol);
    f_vol = strrep(f_vol, '.nii.gz', '.nii');
    F_cleanup = {f_vol};
end

[~, h_vol] = hb_nii_load(f_vol, 'JustGetHeader', 1);
N_vol = length(h_vol);
if not(isempty(opts.FramesToMap))
    assert(length(opts.FramesToMap)<=N_vol);
    assert(all(ismember(opts.FramesToMap, 1:N_vol)));
    error('code');
end

if ischar(hemi)
    assert(ismember(hemi, {'lh','rh'}));
    hemis = {hemi};
elseif iscell(hemi)
    hemis = hemi;
end
Nh = length(hemis);

assert(Nh<=2);

f_surfs = struct;
osn = opts.OutputSurfaceName;
if isempty(osn)
    [p_vol,n_vol,e] = fileparts(f_vol);
    assert(isequal(e, '.nii'));
    for iH=1:Nh
        hemi = hemis{iH};
        n = sprintf('%s.%s.%s.gii', hemi, n_vol, opts.SurfaceSpace);
        f_surfs.(hemi) = fullfile(p_vol, n);
    end
else
    assert(not(isempty(osn)));
    assert(or(ischar(osn),isstruct(osn)), 'unrecognized input format');
    if Nh==1
        if ischar(osn)
            f_surfs.(hemis{1}) = osn;
        else
            assert(isfield(osn, hemis{1}));
            f_surfs = osn;
        end
    else
        assert(isstruct(osn));
        assert(isfield(osn, hemis{1}));
        assert(isfield(osn, hemis{2}));
        f_surfs = osn;
    end
end

switch opts.SurfaceSpace
    case {'fsaverage', 'fsLR'}
        % duplicate fsaverage_hb to opts.dir_subjs.
        % then cleaned up later below.

        d_fsavg_tmp = fullfile(fsopts.dir_subjs, 'fsaverage_hb');
        copyfile(fsopts.dir_fsavg_hb, d_fsavg_tmp);

        if isempty(F_cleanup)
            F_cleanup = {d_fsavg_tmp};
        else
            F_cleanup(length(F_cleanup)+1) = d_fsavg_tmp;
        end
end
end

%==========================================================================
function docleanup(F)
N = length(F);
for k=1:N
    f = F{k};
    if isfolder(f)
        rmdir(f, 's');
    else
        delete(f);
    end
end
end

%==========================================================================
function runcmd(cmd,msg,f)
[sts,log] = system(cmd);
if ~isempty(f)
    if ~exist(f, 'file')
        fprintf('\n.File not generated: %s', f);
        fprintf('\n.log:');
        log %#ok<NOPRT> 
    end
end
if sts~=0
    sprintf('*** %s log ***\n\n', msg);
    log %#ok<NOPRT> 
    error('%s error.',msg);
end
end

