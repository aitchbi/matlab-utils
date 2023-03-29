function f_o = hb_reslice_vol(f_i,f_r,interp,f_o,silentMode)
%HB_RESLICE_VOL Coregisteres an input volume to a given reference volume,
% and then reslices the volume so that the volume dimentions match,
% resulting in a one-to-one correspondence between the voxels of the output
% volume and the reference volume. The new resampled volume is written to
% the directory of the input volume, unless name of output file specified.
%
% Inputs:
%   f_i: file to reslice; full path.
%   f_r: file to use as ref for resolution & coregistration; full path.
%   interp: interpolation order; see spm_reslice.m for details.
%   f_o: name of output file name to save; full path.
%   silentMode: skip outputing info & also block SPM banners.
%
%   NOTE: if f_i is on a read-only drive, f_o has to be specified. 
%
% Outputs:
%   f_o: resliced file; full path. 
% 
%
% Example usage:
%   h_o = hb_reslice_vol(f_i,f_r);
%   hb_reslice_vol(f_i,f_r,1,f_o);
%
% Dependencies:
%   hb_gunzip.m
%   spm_run_coreg_hb.m
%   spm_reslice_hb.m
%   spm_run_coreg.m (SPM12)
%   spm_reslice.m (SPM12)
%
% See also:
%   hb_resample_vol.m
%
% Hamid Behjat

if ~exist('interp','var') || isempty(interp)
    interp = 0;
end
if all(interp~=[0,1])
    error('Interpolation must be either 0 or 1')
end
if ~exist('f_o','var')
    f_o = [];
end
if ~exist('silentMode','var') || isempty(silentMode)
    silentMode = false;
end

if ~silentMode
    fprintf('\n Interpolation order used for reslicing: %d \n',interp); 
end

if isstruct(f_i)
    f_i = f_i.fname;
end

if ~exist(f_r,'file') 
    hb_gunzip(f_r);
    cleanup_r = true;
else
    cleanup_r = false;
end

if isempty(f_o)
    h_r = spm_vol(f_r);
    d = abs(diag(h_r(1).mat));
    if ~isequal(d(1),d(2),d(3)); error('fishy.'); end
    res_o = d(1);
    d = sprintf('%04d',res_o*1e3);
    [p_i,n_i] = fileparts(f_i);
    f_o = fullfile(p_i,[n_i,'_res',d,'.nii']);
end

if isequal(fileparts(f_i),fileparts(f_o))
    diff_io = false;
else
    diff_io = true;
end

if ~exist(f_i,'file') 
    hb_gunzip(f_i);
    cleanup_i = true;
else
    cleanup_i = false;
end

job.ref = {strcat(f_r,',1')};
if length(spm_vol(f_i))==1 
    job.source = {strcat(f_i,',1')};
else 
    job.source = {f_i};
end
job.roptions.interp = interp;
job.roptions.wrap = [0 0 0];
job.roptions.mask = 0;
job.roptions.prefix = 'tmp_';

if diff_io
    job.roptions.writedirectory = fileparts(f_o);
    if silentMode
        out = spm_run_coreg_hb(job,true);
    else
        out = spm_run_coreg_hb(job);
    end
else
    if silentMode
        out = spm_run_coreg_hb(job,true);
    else
        out = spm_run_coreg(job);
    end
end

f_tmp = strsplit(out.rfiles{1},',');
movefile(f_tmp{1},f_o); % rename tmp_*.nii to name{f_o}

if cleanup_i
    delete(f_i);
end

if cleanup_r
    delete(f_r);
end
