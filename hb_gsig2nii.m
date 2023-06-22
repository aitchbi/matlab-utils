function hb_gsig2nii(S, G, f_o, varargin)
% HB_GSIG2NII write graph signals associated to a voxel-wise graph as a
% nifti file(s).
%
% Inputs:
%   S: graph signal matrix, one signal per column.
%
%   G: (graph) structure, with the fields "indices", "dim", and "mat"
%   (optional); "indices" gives the indices in the 3D volume associated to
%   the graph, one index per graph vertex; "dim" (1 x 3 array) is the
%   dimension of the volume in which the graph was created (should match
%   that of ref), and "mat" (4 x 4 matrix) is an affine that specifies
%   voxel dimensions, and also how voxel coordinates can be transformed to
%   world cordinates. "dim" and "mat" should match that of ref, but only
%   "dim" will be stricly verified, and only a warninig will be issued if
%   field "mat" is misisng, and thus cannot be verified.
%
%   f_o: absolute path of nifti file to be written.
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
%   ReferenceFile: either (a) the absolute path of nifti file or (b) a
%   structure with 3 fileds: dim, mat, dt; see NOTE1 at the end of this
%   function. The specified dimentions and data format will be used as
%   reference to write output file. If ReferenceFile specified, fileds
%   "dim", "mat" and "dt" will be used, the first to to verify match
%   between the graph and reference file, and the third one to specify the
%   data type for writing out the file.
%
%   NanValuesOutsideMask: logical (true or false) By default, values
%   outside the mask (voxels not associated with graph vertices) are assign
%   the value NaN, which enables i) better visulaisation (e.g. in FSLeyes)
%   and ii) verification of non-graph-signal voxels. To assign the value 0
%   instead of NaN set this input to false. [default: true] 
%
%   DataType: one the following values: 2, 4, 8, 16, 64, 256, 512, or 768.
%   See NOTE1 at the end of this function. Default: 64.
%
% NOTE: If neither ReferenceFile nor DataType specified, datatype 64 is
% used by deault. If only ReferenceFile specified, its datafile will be
% used. If only DataType specified, this datafile will be used. If both
% ReferenceFile and DataType specfied, the value specidfed by DataType will
% be used.
%
% Examples:
% [1] basic call, used dim & mat from G + default data type:
%     hb_gsig2nii(S, G, f_o);
%
% [2] specifying the data type to use when writing the output file:
%     hb_gsig2nii(__, 'DataType', 16);
%
% [3] a nifti file per graph signal rather than a 4D nifti:
%     hb_gsig2nii(__, 'WriteSeperateFileForEachSignal', true);
%
% [4] assign 0 to values outside mask instead of NaN:
%     hb_gsig2nii(__, 'NanValuesOutsideMask', false);
%
% [5] specifying a reference file as a file:
%     f_r = '/Users/x/y/z.nii.gz';
%     hb_gsig2nii(__, 'ReferenceFile', f_r);
%
% [6] specifying a reference file as a structure:
%     h_r = struct('dim', [100 100 100], 'mat', eye(4), 'dt', [16 0]);
%     hb_gsig2nii(__, 'ReferenceFile', h_r);
%
% Dependencies: 
% . SPM12 package
%
% Hamid Behjat

% datatype 64 used by default
d = inputParser;
addParameter(d,'WriteSeperateFileForEachSignal', false); 
addParameter(d,'ReferenceFile', []); 
addParameter(d,'NanValuesOutsideMask', true); 
addParameter(d,'DataType', []); 
addParameter(d,'Verbose', true);
parse(d,varargin{:});
opts = d.Results;

%-Check inputs. 
[f_o, f_ogz] = chkinputs(S, G, f_o);

%-Reference header.
[DataDim, DataMat, DataType] = getdmd(G, opts);

%-Prepare output file header.
h_o = gethout(DataDim, DataMat, DataType, f_o, opts);

%-Write volume(s).
writevols(S, G, h_o, opts)

%-Cleanup. 
cleanup(f_o, f_ogz);

end





%==========================================================================
function h_o = gethout(DataDim, DataMat, DataType, f_o, opts)
h_o = struct;
h_o.dim = DataDim;
h_o.mat = DataMat;
h_o.dt  = DataType;
if not(opts.WriteSeperateFileForEachSignal)
    h_o.fname = f_o;
    h_o = spm_create_vol(h_o);
end
end


%==========================================================================
function [f_o, f_ogz] = chkinputs(S, G, f_o)
assert(length(G.indices)==size(S,1));
if endsWith(f_o, '.gz')
    f_ogz = f_o;
    f_o = strrep(f_o, '.gz', '');
else
    f_ogz = [];
end
assert(endsWith(f_o,'.nii'), 'f_o: absolute path of nifti file expected.');
[~,~] = mkdir(fileparts(f_o));
end

%==========================================================================
function [DataDim, DataMat, DataType] = getdmd(G, opts)
% dmd: dim, mat, datatype
if isempty(opts.ReferenceFile)
    DataDim = G.dim;
    DataMat = G.mat;
    DataTypeRefFile = [];
else
    ref = opts.ReferenceFile;
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
    
    % verify dimentions
    msg = 'Inputs G & ref are not matched.';
    assert(isequal(G.dim, h_r.dim), msg);
    
    % also verift affines, if given in G
    if isfield(G,'mat')
        assert(all(abs(G.mat-h_r.mat)<1e-6, 'all'), msg);
    else
        d1 = 'Graph structure missing filed ".mat"; ';
        d2 = 'match between G and ref file not verified.';
        msg = [d1 d2];
        warninig(msg);
    end
    DataDim = h_r.dim;
    DataMat = h_r.mat;
    DataTypeRefFile = h_r.dt(1);
end

% . If both opts.DataType and opts.ReferenceFile given, the former is used. 
% . If neither opts.DataType nor opts.ReferenceFile given 64 is used. 
if isempty(opts.DataType)
    if isempty(DataTypeRefFile)
        DataType = [64 0];
    else
        DataType = [DataTypeRefFile 0];
    end
else
    DataType = [opts.DataType 0];
end
end

%==========================================================================
function cleanup(f_o, f_ogz)
if not(isempty(f_ogz))
    gzip(f_o);
    delete(f_o);
end
end

%==========================================================================
function writevols(S, G, h_o, opts)
Ns = size(S,2);
for k=1:Ns
    if opts.NanValuesOutsideMask
        v_o = nan(h_o.dim);
    else
        v_o = zeros(h_o.dim);
    end
    v_o(G.indices) = S(:,k);
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