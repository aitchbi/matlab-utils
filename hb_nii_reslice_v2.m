function [f_o, M] = hb_nii_reslice_v2(f_i, f_r, f_o, varargin)
%HB_NII_RESLICE_V2 coregisteres an input volume to a given reference volume,
% and then reslices the volume so that the volume dimentions match,
% resulting in a one-to-one correspondence between the voxels of the output
% volume and the reference volume. the new resampled volume is written to
% the directory of the input volume, unless name of output file specified.
%
% inputs:
%   f_i: file to reslice; full path. see NOTE 1. 
%
%   f_r: file to use as ref for resolution & coregistration; full path.
%
%   f_o: output file name to be saved; full path.
%
%   name-value pair arguments:
%   InterpolationOrder: interpolation order; see spm_reslice.m for details.
%
%   Verbose: to skip outputing info + block SPM banners; default: false.
%
%   InputFilesInReadOnlyDir: 1x2 logical array, specifying
%   wherther f_i and/or of f_r are in read-only directories, respectively;
%   default: [false false].
%
%   RegisterThenReslice: a logical; if true, f_i will first be registered
%   to f_r, and then, it will be resliced; default: false.
%
%   RegisterationMatrix: precomputed registration matrix, as computed
%   and output from spm_run_coreg_hb.m; note: if RegisterationMatrix is
%   input, there is no need to set RegisterThenReslice as that is implied
%   to be true.
%
% outputs:
%   f_o: resliced file; full path. 
%
%   M: estimated registration matrix; [] if RegisterThenReslice=false
%
%   NOTE 1
%   to reslice multiple files, input f_i as a cell array of file
%   paths; the first file will be used for regiteration to f_r (e.g. an
%   anatomical) whereas the other files will be resliced based on the
%   header-of/registeration-obtained-based-on the first image. 
% 
%   this option, i.e., multiple input files, makes most sense to use when
%   RegisterThenReslice = true, in which case, a registeration will be
%   obtained based on the first image, and will then be used to reslice all
%   the files (multiple anatomical files, functional, PET, etc.; whatever
%   is originally in register with the first image).
%  
%   alternatively, a precomputed registration matrix can be input
%   (RegisterationMatrix) and if so, no registration will be estimated,
%   instead, the precompute matrix will be used.
%
% example useage: 
%
% [Exp-1] if both input files (f_i and f_r) are in writable directory: 
%
%         hb_nii_reslice(f_i, f_r, f_o);
%
% [Exp-2] if eg f_i in writable directory but f_r in read-only directory: 
%          
%         hb_nii_reslice(f_i, f_r, f_o, 'InputFilesInReadOnlyDir', [false true]);
%
% [Exp-3] first register, then reslice:
%         hb_nii_reslice(f_i, f_r, f_o, 'RegisterThenReslice', true); 
% 
% [Exp-4] use precomputed registratioin matrix for reslicing. in eg belwo
% it is assumed that f_i1 and f_i2 are in register, so the registration
% estimated for taking f_i1 to f_r can be used for mapping f_i2 to f_r;
% this is better approach than recomputing the registration for mapping
% f_i2 to f_r, not only computational-wise but more importantly to reduce
% bias, ensuring f_i1 and f_i2 remain in exact register when mapped to the
% new space.
%  
% [~,~,M] = hb_nii_reslice(f_i1, f_r, f_o1, 'RegisterThenReslice', true); 
% hb_nii_reslice(f_i2, f_r, f_o2, 'RegisterationMatrix', M); 
%
% dependencies:
%   SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12
%   spm_run_coreg_hb.m [*]
%   spm_reslice_hb.m [*]
%   [*] https://github.com/aitchbi/matlab-utils/tree/main/spm_modified
%
% h behjat

p = inputParser;
addParameter(p,'InterpolationOrder', 1);
addParameter(p,'Verbose', false);
addParameter(p,'InputFilesInReadOnlyDir', [false false]); % [f_i f_r]
addParameter(p,'RegisterThenReslice', []); % don't set to false here; see below
addParameter(p,'RegisterationMatrix', []); % n/a if RegisterThenReslice specified by used =false
parse(p,varargin{:});
opts = p.Results;

assert(ismember(opts.InterpolationOrder, [0 1]), 'supported interpolation orders: 0 or 1');

if opts.Verbose
    fprintf('\n interpolation order used for reslicing: %d \n', opts.InterpolationOrder); 
end

if ~isempty(opts.RegisterationMatrix)
    assert(ismatrix(opts.RegisterationMatrix));
    assert(isequal(size(opts.RegisterationMatrix), [4 4]));
    if isempty(opts.RegisterThenReslice)
        opts.RegisterThenReslice = true;
    else
        d1 = 'fishy; registration matrix is input';
        d2 = ', whereas RegisterThenReslice is input as false';
        msg = [d1, d2];
        assert(opts.RegisterThenReslice, msg);
    end
end

assert(or(ischar(f_i), iscell(f_i)));
% char: single file path
% cell: array of file paths

if iscell(f_i)
    if length(f_i)==1
        f_i = f_i{1};
        f_others = [];
        N_others = 0;
    else
        N_others = length(f_i)-1;
        f_others = f_i(2:end);
        f_i      = f_i{1};
    end
else
    N_others = 0;
    f_others = [];
end

if N_others>0
    assert(iscell(f_o));
    assert(length(f_o)==N_others+1);
    f_others_o = f_o(2:end);
    f_o = f_o{1};
end

if isequal(fileparts(f_i), fileparts(f_o))
    diff_io = false;
else
    diff_io = true;
end

if endsWith(f_o, '.nii')
    
    GzipOutput = false;
    
elseif endsWith(f_o, '.gz')
    
    f_o = strrep(f_o, '.gz', '');

    GzipOutput = true;
end

%-use temporary working directory & handle gzips.
%--------------------------------------------------------------------------
if any(opts.InputFilesInReadOnlyDir)
    
    TWD = fullfile(...
        fileparts(f_o),...
        sprintf('hb_nii_displace_%s', get_randtag));
    
    [~,~] = mkdir(TWD);
    
    if opts.InputFilesInReadOnlyDir(1)
        f_i = cp2twd(f_i, TWD);
    end
    
    if opts.InputFilesInReadOnlyDir(2)
        f_r = cp2twd(f_r, TWD);
    end
else
    TWD = [];
end

[f_r, ~, cleanup_r] = handlegzip(f_r);

[f_i, ~, cleanup_i] = handlegzip(f_i);

%-define & run job. 
%--------------------------------------------------------------------------
job.ref = {strcat(f_r,',1')};

if length(spm_vol(f_i))==1 
    job.source = {strcat(f_i,',1')};
else
    job.source = {f_i};
end

job.roptions.interp = opts.InterpolationOrder;
job.roptions.wrap = [0 0 0];
job.roptions.mask = 0;
job.roptions.prefix = sprintf('%s_',get_randtag);

if diff_io
    job.roptions.writedirectory = fileparts(f_o);
end

if opts.RegisterThenReslice
    job.eoptions.cost_fun = 'nmi'; % *
    % *: this is the default cost func; just specified as a flag to have
    % spm_run_coreg performs registration before reslicing.
else
    assert(isempty(opts.RegisterationMatrix), ...
        'fishy; reg matrix is input whereas no registration is requested');
end

if N_others>0
    job.other = f_others;
end

[out, M] = spm_run_coreg_hb(job, not(opts.Verbose), opts.RegisterationMatrix);

f_tmp = strsplit(out.rfiles{1}, ',');

movefile(f_tmp{1}, f_o); % rename tmp_*.nii to name{f_o}

if N_others>0
    for k=1:N_others
        f_tmp = strsplit(out.rfiles{1+k}, ',');
        %if N_others==1
        %    movefile(f_tmp{1}, f_others_o);
        %else
            movefile(f_tmp{1}, f_others_o{k});
        %end
    end
end

%-cleanups.
%--------------------------------------------------------------------------
if cleanup_i
    delete(f_i);
end

if cleanup_r
    delete(f_r);
end

if GzipOutput
    gzip(f_o);
    delete(f_o);
    f_o = [f_o, '.gz'];
end

if ~isempty(TWD)
    rmdir(TWD,'s');
end
end

%==========================================================================
function [f, f_gz, cleanup] = handlegzip(f)
if endsWith(f,'.gz')
    f_gz = f;
    f = strrep(f,'.gz','');
else
    f_gz = [f,'.gz'];
end
if exist(f,'file')
    f_gz = [];
    cleanup = false;
else
    assert(...
        logical(exist(f_gz,'file')),...
        sprintf('GZIP file does not exist: %s',f_gz));
    gunzip(f_gz);
    cleanup = true;
end
end

%==========================================================================
function f2 = cp2twd(f1,TWD)
% copy file to TWD if in read-only dir. 
d = strrep(f1, fileparts(f1), TWD);
while 1 % NOTE
    f2 = strrep(d, '.nii', ['_', get_randtag, '.nii']);
    if ~exist(f2,'file')
        break
    end
end
sts = copyfile(f1,f2);
assert(sts==1);
% NOTE: random tag for robustness, e.g. in parallel runs.
end

%==========================================================================
function t = get_randtag
t = sprintf('tmp%d',round(rand*1e12));
end
