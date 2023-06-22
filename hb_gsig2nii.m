function hb_gsig2nii(S, I, f_o, ref, varargin)
% 
% Inputs:
%
%   S: graph signal matrix, one signal per column.
%
%   I: indices of voxels in 3D volume corresponding to graph vertices
%
%   f_o: absolute path of nifti file to be written.
%
%   ref: either (a) the absolute path of nifti file or (b) a
%   structure with 3 fileds: dim, mat, dt; see NOTE1 at the end of this
%   file. The specified dimentions and data format will be used as
%   reference to write output file.
%
%   Format for "f_o": '/x/y/z.e' where '.e' is '.nii' or '.nii.gz'. Same for
%   "ref" if given as a filepath not a structure.  
%
%   Name-Value Arguments:
%   WriteSeperateFileForEachSignal: given multiple graph signals, write all
%   as a single 4D nifti file, or each seperately, i.e., a set of 3D
%   niftis. For the latter, the given name will be adapted to reflect the
%   order of columns in S, i.e, e.g. first signal (first column in S) will
%   be tagged as *_0001.nii.
%
% Examples:
% [1] f_r = '/Users/x/y/z.nii.gz';
%     hb_gsig2nii(S, I, f_o, f_r);
%
% [2] h_r = struct('dim', [100 100 100], 'mat', eye(4), 'dt', [16 0]);
%     hb_gsig2nii(S, I, f_o, h_r);
%
% Dependencies: 
% . SPM12 package
%
% HB

d = inputParser;
addParameter(d,'WriteSeperateFileForEachSignal', false); 
addParameter(d,'Verbose', true);
parse(d,varargin{:});
opts = d.Results;

if endsWith(f_o, '.gz')
    f_ogz = f_o;
    f_o = strrep(f_o, '.gz', '');
else
    f_ogz = [];
end
assert(endsWith(f_o,'.nii'), 'f_o: absolute path of nifti file expected.');

assert(length(I)==size(S,1));

Ns = size(S,2);

[~,~] = mkdir(fileparts(f_o));

%-Reference header.
if ischar(ref)
    f_r = ref;
    if not(exist(f_r,'file'))
        f_rgz = [f_r, '.gz'];
        if exist(f_rgz,'file')
            f_r = f_rgz;
        end
    end
    h_r = spm_vol([f_r, ',1']);
else
    assert(isstruct(ref));
    assert(isfield(ref, 'dim'));
    assert(isfield(ref, 'mat'));
    assert(isfield(ref, 'dt'));
    h_r = ref;
end

%-Prepare header.
h_o = struct;
h_o.dim = h_r.dim;
h_o.mat = h_r.mat;
h_o.dt  = h_r.dt;

if not(opts.WriteSeperateFileForEachSignal)
    h_o.fname = f_o;
    h_o = spm_create_vol(h_o);
end

% write volumes
for k=1:Ns
    v_o = zeros(h_o.dim);
    v_o(I) = S(:,k);
    if opts.WriteSeperateFileForEachSignal
        h_o.fname = strrep(f_o, '.nii', sprintf('_%04d.nii', k));
    else
        h_o.n(1) = k;
    end
    spm_write_vol(h_o, v_o);
    if opts.Verbose && Ns>1
        showprgs(k, Ns);
    end
end

if not(isempty(f_ogz))
    gzip(f_o);
    delete(f_o);
end

end

%==========================================================================
function showprgs(n,N)
l = numel(num2str(N));
if n==1
    fprintf('\n..Mapping volume using displacement map.. ');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end

%-NOTES--------------------------------------------------------------------
%
%-NOTE1.
% Structure with either (or all) of these three fields: 
% .dim: 1x3 positive integer array 
% .mat: 4x4 numeric array
% .dt: 1x2 array, where second element is typically 0, and the first value
% specifies the data type*, which will in turn notably determine the size of
% the file that is saved, i.e., number fo bits used to store each voxel
% value.
% 
% * any of the values from "types" below, e.g. 4 or 64.
%
%-Extract from spm_type.m--------------------------------------------------
%prec   = {'uint8','int16','int32','float32','float64','int8','uint16','uint32'};
%conv   = {@uint8,@int16,@int32,@single,@double,@int8,@uint16,@uint32};
%types  = [    2      4      8   16   64   256    512    768];
%maxval = [2^8-1 2^15-1 2^31-1  Inf  Inf 2^7-1 2^16-1 2^32-1];
%minval = [    0  -2^15  -2^31 -Inf -Inf  -2^7      0      0];
%nanrep = [    0      0      0    1    1     0      0      0];
%bits   = [    8     16     32   32   64     8     16     32];
%--------------------------------------------------------------------------