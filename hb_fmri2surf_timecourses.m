function S = hb_fmri2surf_timecourses(f_fmri, FS, varargin)
% HB_FMRI2SURF take a 4D fMRI volume to FreeSurfer (FS) space and extract
% time-courses for surface vertices.
%
% inputs:
%   f_rs: rs_processed_censored_32bit.nii
%
%     FS: structure with fields 'f_ribbon', 'surfspace', 'fsopts', e.g.: 
%
%       FS = struct;
%       FS.f_ribbon              = fullfile(param.hcp.dnames.t1w, param.subject, 'mri', 'ribbon.nii');
%       FS.surfspace             = 'subject';
%       FS.fsopts.ID             = param.subject;
%       FS.fsopts.dir_subjs      = param.hcp.dnames.t1w;
%       FS.fsopts.sh_vol2surf    = fullfile(param.d_hb, 'hb_vol2surf.sh');
%       FS.fsopts.dir_fsavg_hb   = param.dir_fsaverage_hb;
%       FS.fsopts.dir_freesurfer = param.paths.FreeSurfer;
%
%  f_out: output file; rs file match to fs-space parcellation.
%
% output:
%   S: structure with fields.. see below. 
%
% dependencies:
%   .SPM12 toolbox
%   .hb_reslice_vol.m
%
% h behjat

d = inputParser;
addParameter(d,'ResliceToMatchFreeSurferRibbon', false); 
parse(d,varargin{:});
opts = d.Results;

h_fmri   = spm_vol(f_fmri);
dim_fmri = h_fmri(1).dim;
mat_fmri = h_fmri(1).mat;
dt_fmri  = h_fmri(1).dt;
N_t      = length(h_fmri);
tag = get_randtag;

d = strrep(f_fmri, '.gz', '');
f_tmp1 = strrep(d, '.nii', sprintf('_%s_step1.nii',tag));
f_tmp2 = strrep(d, '.nii', sprintf('_%s_step2.nii',tag));

[p,n,e] = fileparts(f_tmp2);
assert(isequal(e,'.nii'));
f_gii = struct;
f_gii.lh = fullfile(p, ['lh.',n,'.gii']);
f_gii.rh = fullfile(p, ['rh.',n,'.gii']);

h0       = struct;
h0.fname = f_tmp1;
h0.dim   = dim_fmri;
h0.mat   = mat_fmri;
h0.dt    = dt_fmri;
h0 = spm_create_vol(h0);

%-output structure.
%--------------------------------------------------------------------------
S = struct;
S.fmri.fname       = f_fmri;
S.fmri.header.dim  = dim_fmri;
S.fmri.header.mat  = mat_fmri;
S.fmri.data.lh     = [];
S.fmri.data.rh     = [];
S.fmri.data.ctx    = [];
S.fmri.data.subctx = [];

S.fmri.data.format = 'number-of-vertices x number-of-frames';
S.fmri.data.descrip = {
    sprintf('fmri.fname was: %s', f_fmri)
    '1) projected to surface (via hb_vol2surf)'
    '2) rsfmri vertex values in lh/rh urface extracted'
    '3) each vector (frame/time-point) saved as a column in matrix.'
    };
S.surface.space = FS.surfspace;
S.surface.N_vertices_lh = [];
S.surface.N_vertices_rh = [];

%-run.
%--------------------------------------------------------------------------

for iT=1:N_t
    
    assert(isequal(h_fmri(iT).dim, dim_fmri));
    
    assert(isequal(h_fmri(iT).mat, mat_fmri));
    
    assert(isequal(h_fmri(iT).dt, dt_fmri));
    
    %-extract single frame-------------------------------------------------
    v = spm_read_vols(h_fmri(iT));
    
    spm_write_vol(h0,v);
    
    %-reslice to match FS ribbon if needed---------------------------------
    if opts.ResliceToMatchFreeSurferRibbon

        RegisterThenReslice = false;

        hb_nii_reslice(f_tmp1, FS.f_ribbon, 1, f_tmp2, true, [], RegisterThenReslice);

        if iT==1
            sts = verifreg(f_tmp2, FS.f_ribbon, 'ribbon', false);
            if sts==0
                fprintf('\n.Trying again. First register, then reslice.')
                delete(f_tmp2);
                RegisterThenReslice = true;
                hb_nii_reslice(f_tmp1, FS.f_ribbon, 1, f_tmp2, true, [], RegisterThenReslice);
                verifreg(f_tmp2,FS.f_ribbon, 'ribbon', true);
                fprintf('\n.Register then reslice worked!');
                fprintf('\n.Note: upto 10s/frame longer process time.');
            end
        end

    else

        f_tmp2 = f_tmp1;
        % for HCP, and for the version of fMRI data i work with (i.e.,
        % those brought from MNI spcae to T1w/Diffusion space), i think the
        % below registration and reslicing is not needed; registration for
        % sure should not be needed since the input f_fmri (either being
        % defined in T1w space of Diffusion space) are both in register;
        % but reslicing is needed only if the projection to surface
        % pipeline really requires the voume file to be in voxel-to-voxel
        % register with the ribbon file, which i don't think it is true and
        % it does not make sense to be so as this is someting that can (and
        % should/could be) handled in the freesurfer function; but if this
        % is really needed, then we should reslice the input voume to match
        % voxel-by-voxel to the resolution of the FS ribbon, that is 1mm3
        % voxels and quite a specfic bounding box.
    end

    %-project to surface---------------------------------------------------
    
    hb_vol2surf(f_tmp2, {'lh','rh'}, FS.fsopts, ...
        'OutputSurfaceName', f_gii, ...
        'SurfaceSpace', FS.surfspace);
    
    assert(exist(f_gii.lh, 'file'), 'f_gii.lh missing');

    assert(exist(f_gii.rh, 'file'), 'f_gii.rh missing');

    %-extract urface vertex values-----------------------------------------
    [maplh, maprh] = loadgii(f_gii, iT);

    if iT==1
        Nlh = length(maplh);
        Nrh = length(maprh);
        S.surface.N_vertices_lh = Nlh;
        S.surface.N_vertices_rh = Nrh;
        S.fmri.data.lh  = zeros(Nlh, N_t);
        S.fmri.data.rh  = zeros(Nrh, N_t);
    end

    S.fmri.data.lh(:,iT) = maplh;
    S.fmri.data.rh(:,iT) = maprh;

    showprgs(iT,N_t);
end

delete(f_tmp1);
delete(f_tmp2);

delete(f_gii.lh);
delete(f_gii.rh);
end

%==========================================================================
function t = get_randtag
t = sprintf('tmp%d',round(rand*1e16));
end

%==========================================================================
function runcmd(cmd,iT)
[sts, log] = system(cmd);

if sts==0
    if iT==1
        fprintf('\n..ANTs successful on first frame.\n');
        log
    end
else
    sprintf('*** ANTS''s log ***');
    log %#ok<*NOPRT>
    error('ANTS error.');
end
end

%==========================================================================
function showprgs(n,N)
l = numel(num2str(N));
if n==1
    fprintf('\n..projecting fMRI vol to surf.. ');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end

%==========================================================================
function [maplh, maprh] = loadgii(f_gii, iT)
d = gifti(f_gii.lh);
maplh = d.cdata;
d = gifti(f_gii.rh);
maprh = d.cdata;
if iT==1
    assert(size(maplh,1)>size(maplh,2));
    assert(size(maplh,2)==1);
    assert(size(maprh,1)>size(maprh,2));
    assert(size(maprh,2)==1);
end
end

%==========================================================================
function sts = verifreg(f_fmri,f_ab,type,ThrowError)
th = 0.05; % percenatge of mismatch, fraction of 1  
[v_rs,h_rs] = hb_nii_load(f_fmri);
[v_ab,h_ab] = hb_nii_load(f_ab);
hb_nii_verify_space_match(h_rs, h_ab);
v_rs(isnan(v_rs)) = 0;
v_ab(isnan(v_ab)) = 0;
v1 = imfill(logical(v_rs), 'holes');
v2 = imfill(logical(v_ab), 'holes');
v3 = and(v1,v2);
n1 = nnz(v1);
n2 = nnz(v2);
n3 = nnz(v3);
ndiff1 = n1-n3;
ndiff2 = n2-n3;
assert(ndiff1>=0, 'fishy');
assert(ndiff2>=0, 'fishy');
p1d = ndiff1/n1;
p2d = ndiff2/n2;
t = 'percentage of';
disp(' ');
fprintf('----------registeration check.');
fprintf('\n.registeration Verification Report:');
fprintf('\n..RS file: %s', f_fmri);
fprintf('\n..anatomical file: %s', f_ab);
switch type
    case 't1'
        fprintf('\n..%s RS voxels outside anatomical: %0.1f%% (P1)',t,p1d*100);
        fprintf('\n..%s anatomical voxels outside RS: %0.1f%% (P2)',t,p2d*100);
        fprintf('\n..pass criteria: P1, P2 or both should be less than %d%%',th*100);
        if or(p1d<=th, p2d<=th)
            fprintf('\n..registeration between RS & anatomical image is OK.');
        else
            error('registeration between RS & anatomical image is fishy.');
        end
    case 'ribbon'
        fprintf('\n..%s ribbon voxels without RS value: %0.1f%% (P)',t,p2d*100);
        fprintf('\n..pass criteria: P should be less than %d%%',th*100);
        if p2d<=th
            fprintf('\n..registeration between RS & ribbon is OK.');
            sts = 1;
        else
            sts = 0;
            if ThrowError
                error('registeration between RS & ribbon is fishy.');
            else
                fprintf('\n..registeration between RS & ribbon is fishy.');
            end
        end
end
fprintf('\n----------registeration check.');
disp(' ');
end
