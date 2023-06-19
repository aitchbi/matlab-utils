function [S, f_reg] = hb_nii2gsig(f, G, varargin)
% HB_NII2GSIG extract graph signal from given input nifti file.
%
% Inputs:
%   f: full path of nifti file
%   G: graph structure with fields indices etc.  
%
%   Name-Value Pair Arguments:
%   GraphRefNii: if fields .dim and .mat not given in G, specify a
%   reference file which is in register with the space in which the graph
%   was defined (T1wSPace or DiffusionSpace depending on the graph);
%   typically the G.f.mask/source file. 
%
%   NiiAndGInRegister: logical (true or false; default: false). If you are
%   sure f and G are in register, set to true to bypass verification.
%
%   ResliceNiiIfNotInRegisterWithG: logical (true or false; default:
%   false). If set to true, GraphRefNii should also be set, in which case f
%   will be resliced to match the space in which G was defiend
%   (GraphRefNii). The resliced file (f_reg) is returned, which is best to
%   be be verified.
% 
% Outputs: 
%   S: graph signals, one per column.
%
% Dependencies: 
% . SPM12 package
% . hb_nii_load.m
%
% HB

%-process optional inputs--------------------------------------------------
funcLogi = @(x) assert(islogical(x));
funcPath = @(x) assert(ischar(x) || isempty(x));
d = inputParser;
addParameter('GraphRefNii', [], funcPath);
addParameter('NiiAndGInRegister', false, funcLogi);
addParameter('ResliceNiiIfNotInRegisterWithG', false, funcLogi);
parse(d,varargin{:});
opts = d.Results;
%--------------------------------------------------------------------------

if endsWith(f, '.gz')
    fgz = f;
    f = strrep(f, '.nii.gz', '.nii');
else
    fgz = [f, '.gz'];
end

if exist(f, 'file')
    CleanUp = false;
else
    if exist(fgz, 'file')
        gunzip(fgz);
        CleanUp = true;
    else
        error('file not found');
    end
end

%-Verify that f & G are in register.
if opts.NiiAndGInRegister
    f_reg = [];
else
    h_f = spm_vol([f, ',1']);
    if isempty(opts.GraphRefNii)
        h_G = struct;
        h_G.dim = G.dim;
        h_G.mat = G.mat;
    else
        h_G = spm_vol(hb_gunzip(opts.GraphRefNii));
    end
    chk1 = isequal(h_f.dim, h_G.dim);
    chk2 = all(abs(h_f.mat-h_G.mat)<1e-6,'all');
    if any([chk1, chk2])
        if opts.ResliceNiiIfNotInRegisterWithG
            fprintf('\n..Reslicing input to match G..');
            f_reg = hb_nii_reslice(f, opts.GraphRefNii, 1);
        else
            errmsg = 'File not in register with space in which graph was defined.';
            assert(chk1, errmsg);
            assert(chk2, errmsg);
            f_reg = [];
        end
    end
end

I = G.indices;
[v, h] = hb_nii_load(f, 'IndicesToLoad', I);

Ng = length(I);

Nv = length(h);

S = zeros(Ng, Nv);
for k=1:Nv
    d = v(:,:,:,k);
    S(:, k) = d(I) ;
end

if CleanUp
    delete(f);
end
end