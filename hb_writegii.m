function hb_writegii(X,S,Y,paths,varargin)
% HB_WRITEGII writes scalar map as gifti based on a reference surface file.
%
% Inputs:
%   X: A vector of length equal to number of vertices in the surface file.
%
%   S: FreeSurfer surface file in curv format; e.g., lh.white
%
%   Y: Output gifti file to write; e.g.,lh.<file-name>.gii
%
%   paths: structure with fields 'freesurfer' and 'sh_asc2gii', which are
%   paths to freesurfer and hb_asc2gii.sh, respectively.
%
% Dependencies:
%   hb_asc2gii.sh [*]
% 
% [*]: give execute permision to hb_asc2gii.sh when using thi function for
% first time. On terminal, cd to the directory where the script is saved
% and then: chmod u+x hb_asc2gii.sh
% 
% Hamid Behjat

p = inputParser;
addParameter(p,'Cleanup', 1);
parse(p,varargin{:});
opts = p.Results;

%-Verify inputs.
dochk(X, S, Y);

%-Calling script path.
d = dbstack('-completenames');
d_scr = fileparts(d(end).file);

%-Temp working directory. 
f_out = Y;
if contains(f_out, filesep)
    d_twd = fileparts(f_out);
else
    d_twd = d_scr;
end
d_twd = fullfile(d_twd, tmptag);
mkdir(d_twd);

%-Read surface.
if contains(S, filesep)
    f_srf = S;
else
    f_srf = fullfile(d_scr, S);
end
vc = read_surf(f_srf);

%-Verify match between data & surface.
N = size(X,1);
assert(size(vc,1)==N);

%-Write data as ascii.
f_in = fullfile(d_twd, [tmptag, '.asc']);
fid = fopen(f_in, 'w');
for k=1:N
    fprintf(fid,...
        '%s %.05f %.05f %.05f %.05f\n',...
        gettag(k-1),...
        vc(k,1),...
        vc(k,2),...
        vc(k,3),...
        X(k));
end
fclose(fid);

%-Convert ascii to gifti.
cmd = sprintf(...
    '%s -i %s -s %s -o %s -r %s',...
    paths.sh_asc2gii,...  % hb_asc2gii.sh
    f_in,...              %i .asc
    f_srf,...             %s *h.white
    f_out,...             %o .gii
    paths.freesurfer);    %r /path/to/freesurfer

[sts,log] = system(cmd);    

if sts~=0
    if contains(lower(log), 'permission denied')
        sprintf('.System log:\n');
        disp(log);
        docleanup(d_twd, opts);
        error('Do: ''chmod u+x hb_asc2gii.sh''');
    else
        sprintf('.FreeSurfer log:\n');
        disp(log);
        docleanup(d_twd, opts);
        error('FreeSurfer error.');
    end
end

%-Cleanup.
docleanup(d_twd, opts);
end

%==========================================================================
function dochk(X,S,Y)
assert(size(X,2)==1, 'Multi-map write not supported.');
if contains(S, filesep)
    d = getname(S);
else
    d = S;
end
assert(startsWith(d, {'lh.', 'rh.'}));
assert(endsWith(Y, '.gii'));
end

%==========================================================================
function t = tmptag
t = sprintf('tmp%d',round(rand*1e6));
end

%==========================================================================
function n = getname(x)
n = strrep(x, fileparts(x), '');
n = strrep(n, filesep, '');
end

%==========================================================================
function docleanup(d,opts)
if opts.Cleanup
    rmdir(d, 's');
end
end

%==========================================================================
function t = gettag(k)
if k<1e3
    t = sprintf('%03d', k);
elseif k<1e4
    t = sprintf('%04d', k);
elseif k<1e5
    t = sprintf('%05d', k);
elseif k<1e6
    t = sprintf('%06d', k);
else
    error('fishy');
end
end