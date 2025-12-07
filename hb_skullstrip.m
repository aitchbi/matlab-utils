function hb_skullstrip(f_t1, f_t1brain, varargin)
% HB_SKULLSTRIP strips away the skull from a T1 image.
%
% h behjat

p = inputParser;
addParameter(p,'Method','FSL');
addParameter(p,'Path_fsl',[]);
parse(p,varargin{:});
opts = p.Results;

if exist(f_t1brain, 'file')
    return;
end

[f_tmp, d_tmp] = duplicatet1(f_t1);

%-strip.
switch opts.Method
    case 'FSL'
        try
            stripskull(f_tmp, f_t1brain, opts.Path_fsl);
        catch ME
            disp(ME);
            if exist(f_t1brain, 'file')
                fprintf('\n.Fishy error in FSL BET, but brain extracted.');
            else
                % once I got this error:
                % '/usr/local/fsl/bin/bet: line 399: 212508 Segmentation fault' [25.09.2023]
                error('Error in FSL BET; brain not extracted.');
            end
        end

    otherwise
        error('extend');
end

%-cleanup.
d = fileparts(f_t1brain);
F = dir(d);
for k=1:length(F)
    n = F(k).name;
    if contains(n, '_tmp_')
        delete(fullfile(d, n));
    end
end
d = strrep(f_t1brain, '.nii', '_mask.nii');
if exist(d, 'file')
    delete(d);
end
delete(f_tmp);
rmdir(d_tmp);
end

%==========================================================================
function stripskull(f_t1, f_t1brain, d_fsl)
if exist(f_t1brain,'file')
    return;
end
if endsWith(f_t1, '.nii.gz')
    gunzip(f_t1);
    f_t1 = strrep(f_t1, '.nii.gz', '.nii');
    cleanupt1nii = true;
else
    cleanupt1nii = false;
end
if endsWith(f_t1brain, '.nii.gz')
    f_t1brain = strrep(f_t1brain, '.nii.gz', '.nii');
end
fprintf('\n..stripping skull.. ');
f_fslbet = hb_fsl_get_func(d_fsl, 'bet');
cmd = sprintf('%s %s %s -B', f_fslbet, f_t1, f_t1brain);
hb_runcmd(cmd,'error in FSL bet skull stripping.');
h = spm_vol(f_t1brain);
v = spm_read_vols(h);
h.descrip = [h.descrip, ' | skull-stripped with FSL bet -B'];
spm_write_vol(h,v);
assert(endsWith(f_t1brain, '.nii'))
gzip(f_t1brain);
delete(f_t1brain);
if cleanupt1nii
    assert(endsWith(f_t1, '.nii'));
    delete(f_t1);
end
end

%==========================================================================
function [f_tmp, d_tmp] = duplicatet1(f)
if endsWith(f, '.gz')
    [p, n, e] = fileparts(f);
    assert(isequal(e, '.gz'));
    assert(endsWith(n, '.nii'));
    n = strrep(n, '.nii', '');
    e = '.nii.gz';
else
    [p, n, e] = fileparts(f);
    assert(isequal(e, '.nii'));
end
d_tmp = fullfile(p, ['hb_skullstrip_', get_randtag]);
f_tmp = fullfile(d_tmp, [n, '_', get_randtag, e]);
mkdir(d_tmp);
copyfile(f, f_tmp);
end

%==========================================================================
function t = get_randtag
t = sprintf('tmp%d',round(rand*1e12));
end

