% [12.06.2026] a stricter (more accurate!) version of
% hb_nii_verify_space_match.m as it appear that there is a need to also
% check the private segment of the nifti header to ensure accurate matching
% of the data in terms of space. in particular, files f1 and f2 may match
% in their h1.mat and h2.mat but may have different h1.private.mat0 and
% h2.private.mat0 and this is why that for some file (eg ADNI prepric rs
% data from Franzmeier) we see that the preproc rs data and the mean rs
% match in their h.mat but when you view them on fesleyes they don't
% overlap and this is because we see a difference between their
% h.private.mat0 which is sort of an additional transformation that is
% looked at when the files are loaded in the viewer.

function [sts, h1, h2, errors] = hb_nii_verify_space_match_v2(f1, f2, varargin)
% HB_NII_VERIFY_SPACE_MATCH verifies match between two nifti files in terms
% of their 3D dimension and coordinate space. that is, whether the two
% vloumes are matched voxel-to-voxel in space. the two inputs can be either
% file paths or nifti headers.
%
% note: if f1, f2, or both are 4D, the assumption is that all frames have a
% matching header, and thus, only the first frame in the AD volume will be
% checked.
% 
% h behjat

d = inputParser;
addParameter(d,'DuplicateThenUnzip', false);
addParameter(d,'ThrowError', false);
parse(d,varargin{:});
opts = d.Results;

h1 = gethead(f1, opts.DuplicateThenUnzip);
h2 = gethead(f2, opts.DuplicateThenUnzip);
chk1 = isequal(h1.dim, h2.dim);
chk2 = all(abs(h1.mat-h2.mat)<1e-2,'all'); % diff <0.01 mm in each element
sts = chk1 && chk2;
if opts.ThrowError
    if not(sts)
        fprintf( ...
            '\n.space mismatch bw files f1 & f2: \n..f1: %s \n..f2: %s', ...
            h1.fname, ...
            h2.fname);
        error('mismatch between files');
    end
end

d1 = h1.private.dat.dim;
d2 = h2.private.dat.dim;
chk3 = isequal(d1(1:3), d2(1:3));
chk4 = all(abs(h1.private.mat - h2.private.mat)<1e-2,'all'); % diff <0.01 mm in each element
chk5 = all(abs(h1.private.mat0 - h2.private.mat0)<1e-2,'all'); % diff <0.01 mm in each element
sts = chk3 && chk4 && chk5;
if opts.ThrowError
    if not(sts)
        fprintf( ...
            '\n.space mismatch bw files f1 & f2 based on info in their header.private : \n..f1: %s \n..f2: %s', ...
            h1.fname, ...
            h2.fname);
        if not(chk3)
            fprintf('\n.mismatch bw dim:');
            fprintf('\n.. h1.private.dat.dim:');
            h1.private.dat.dim(1:3)
            fprintf('\n.. h2.private.dat.dim:');
            h2.private.dat.dim
        end
        if not(chk4)
            fprintf('\n.mismatch bw mat:');
            fprintf('\n.. h1.private.mat:');
            h1.private.mat
            fprintf('\n.. h2.private.mat:');
            h2.private.mat
        end
        if not(chk5)
            fprintf('\n.mismatch bw mat0:');
            fprintf('\n.. h1.private.mat0:');
            h1.private.mat0
            fprintf('\n.. h2.private.mat0:');
            h2.private.mat0
        end
        error('mismatch between files');
    end
end

errors = cell(0,1);
cnt = 0;
if not(chk1)
    cnt = cnt+1;
    errors(cnt) = {'mismatch bw h.dim'};
end
if not(chk2)
    cnt = cnt+1;
    errors(cnt) = {'mismatch bw h.mat'};
end
if not(chk3)
    cnt = cnt+1;
    errors(cnt) = {'mismatch bw h.private.dat.dim'};
end
if not(chk4)
    cnt = cnt+1;
    errors(cnt) = {'mismatch bw h.private.mat'};
end
if not(chk5)
    cnt = cnt+1;
    errors{cnt} = ['mismatch bw h.private.mat0'];
end
end

%==========================================================================
function h = gethead(f,DuplicateThenUnzip)
if ischar(f)
    if 1 % 01.05.2024
        [~,h] = hb_nii_load(f, ...
            'JustGetHeader', true, ...
            'DuplicateThenUnzip', DuplicateThenUnzip);
    else
        if exist(f,'file')
            if endsWith(f, '.nii.gz')
                gunzip(f);
                d = strrep(f, '.gz', '');
                h = spm_vol(d);
                delete(d);
            else endsWith(f, '.nii')
                h = spm_vol(f);
            end
        else
            fgz = [f,'.gz'];
            if exist(fgz, 'file')
                gunzip(fgz);
                h = spm_vol(f);
                delete(f);
            else
                error('Missing file: %s',f);
            end
        end
    end
else
    assert(isstruct(f));
    h = f;
end
if length(h)>1
    h = h(1);
end
end
