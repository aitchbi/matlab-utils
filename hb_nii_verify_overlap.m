function [sts, msg, p1, p2] = hb_nii_verify_overlap(f1,f2,p,type,varargin)
% HB_NII_VERIFY_OVERLAP verifies overlap between image masks e.g. two T1
% images of the same subject after being registered and resliced to match
% voxel-by-voxel. 
%
% inputs:
%   f1: path of file 1.
%
%   f2: path of file 2.
%
%   p: a scalar; overlap percentage in (1 100].
%
%   type: 'both' or 'either' or 'second'; should 'both' images at least
%   have p% overlap, or 'either', or just the 'second'? most stringent
%   check is for 'both', however, that may not be suitable eg if one image
%   is brain masked whereas the other one is not. 
%
% outputs:
%   sts: logicical; verification passed?
%
%   p1: percentage of overlap for f1
%
%   p2: percentage of overlap for f2
%
% dependencies:
% .github/aitchbi/matlab-utils
% .SPM12
%
% examples:
%   sts = hb_nii_verify_overlap(f1, f2, 95, 'either');
%   sts = hb_nii_verify_overlap(f1, f2, 95, 'second');
%
% h behjat

d = inputParser;
addParameter(d, 'MinimumOverlap', []); % n/a if type=='both'
addParameter(d, 'AssumeNansAsZeros', false);
parse(d, varargin{:});
opts = d.Results;

assert(hb_nii_verify_space_match(f1,f2), 'images are not registered');
assert(p>1 && p<=100, 'p should be in (1 100]'); 

v1 = imfill_pvt(f1, opts);
v2 = imfill_pvt(f2, opts);
vo = and(v1,v2);
no = nnz(vo);
p1 = (no/nnz(v1))*100;
p2 = (no/nnz(v2))*100;
c1 = p1>=p;
c2 = p2>=p;

switch type
    case 'both'
        % both should have at least p% overlap
        sts = and(c1,c2);
        if sts
            msg = sprintf('OK. [both have at least %0.1f%% overlap]', p);
        else
            if and(not(c1),not(c2))
                msg = sprintf('not OK. [both files have less than %0.1f%% overlap with each other]', p);
            elseif not(c1)
                msg = sprintf('not OK. [less than %0.1f%% of file 1 overlaps with file 2]', p);
            else
                msg = sprintf('not OK. [less than %0.1f%% of file 2 overlaps with file 1]', p);
            end
        end

    case 'either'
        t1 = sprintf('[at least one has overlap above %0.1f]', p);
        if isempty(opts.MinimumOverlap)
            % at least one should have at least p% overlap
            sts = or(c1,c2);
            if sts
                msg = ['OK. ', t1];
            else
                msg = 'not OK. [neither has overlap above threshold]';
            end
        else
            pm = opts.MinimumOverlap;
            if p1<pm
                sts = false;
                t2 = sprintf('[but file 1 has %0.1f%% overlap with file 2]', p1);
                msg = sprintf('not OK. %s %s', t1, t2);
            elseif p2<pm
                sts = false;
                t2 = sprintf('[but file 2 has %0.1f%% overlap with file 1]', p2);
                msg = sprintf('Not OK. %s %s', t1, t2);
            else
                sts = true;
                t2 = sprintf('[& both above minimum %0.1f%% overlap]', p);
                msg = sprintf('OK. %s %s', t1, t2);
            end
        end

    case 'second'
        % e.g. when file 2 is ribbon and we just want to check that it has
        % good overlap with file 1 that is RS/PET.
        sts = c2;
        if sts
            msg = sprintf('OK. [second has %0.1f%% overlap with first]', p2);
        else
            msg = sprintf('Not OK. [second has %0.1f%% overlap with first]', p2);
        end
end
end

%==========================================================================
function v = imfill_pvt(f, opts)   
try
    v = imfill(logical(hb_nii_load(f)), 'holes');
catch err
    if contains(err.message, 'NaN')
        if opts.AssumeNansAsZeros
            f_tmp = strrep(f, '.nii', '___tmp_hb_nii_verify_overlap.nii');
            f_tmp = strrep(f_tmp, '.gz', '');
            [v, h] = hb_nii_load(f);
            Inan = find(isnan(v));
            v_tmp = v;
            h_tmp = h;
            v_tmp(Inan) = 0;
            h_tmp.fname = f_tmp; 
            spm_write_vol(h_tmp, v_tmp);
            v = imfill(logical(hb_nii_load(f_tmp)), 'holes');
            delete(f_tmp);
        end
    else
        disp(err);
        error('hb_nii_verify_space_match: imfill failed');
    end
end
end
