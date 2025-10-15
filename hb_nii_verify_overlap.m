function [sts, msg, p1, p2] = hb_nii_verify_overlap(f1,f2,p,type,varargin)
% HB_NII_VERIFY_OVERLAP verifies overlap between image masks e.g. two T1
% images of the same subject after being registered and resliced to match
% voxel-by-voxel. 
%
% Inputs:
%   f1: path of file 1.
%
%   f2: path of file 2.
%
%   p: a scalar; overlap percentage in (1 100].
%
%   type: 'both' or 'either' or 'second'; should both images at least have
%   95% overlap or just at least one of them? For a fixed p, 'both' results
%   in a stricter check.
%
% Outputs:
%   sts: logicical; verification passed?
%
%   p1: percentage of overlap for f1
%
%   p2: percentage of overlap for f2
%
% Dependencies:
% .github/aitchbi/matlab-utils
% .SPM12
%
% Example:
% sts = hb_nii_verify_overlap(f1,f2,95,'one');
%
% HB

d = inputParser;
addParameter(d, 'MinimumOverlap', []); % n/a if type=='both'
parse(d, varargin{:});
opts = d.Results;

assert(hb_nii_verify_space_match(f1,f2), 'images are not registered');
assert(p>1 && p<=100, 'p should be in (1 100]'); 
v1 = imfill(logical(hb_nii_load(f1)), 'holes');
v2 = imfill(logical(hb_nii_load(f2)), 'holes');
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
            msg = sprintf('Not OK. [at least one file has less than %0.1f%% overlap]', p);
        end

    case 'either'
        t1 = sprintf('[at least one has overlap above %0.1f]', p);
        if isempty(opts.MinimumOverlap)
            % at least one should have at least p% overlap
            sts = or(c1,c2);
            if sts
                msg = ['OK. ', t1];
            else
                msg = 'Not OK. [neither has overlap above threshold]';
            end
        else
            pm = opts.MinimumOverlap;
            if p1<pm
                sts = false;
                t2 = sprintf('[but file 1 has %0.1f%% overlap with file 2]', p1);
                msg = sprintf('Not OK. %s %s', t1, t2);
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