function [f, f_dlt, sts] = hb_nii_handlegzip(f)
% HB_NII_HANDLEGZIP handles an input nifti file, and returns a version that
% can be loaded, if it exists. If file is in .nii.gz format, or if it is in
% .nii format but only its .nii.gz version exists on disk (in same path a
% input file), the .gz version will be unziped, and the returned file name
% is always .nii (i.e., a file that can be directly loaded).
% 
% Inputs:
%   f: a file path, ending in .nii or .nii.gz. 
% 
% Outputs:
%   f: version of input file that can be loaded, i.e., not in .gz fomrat.
%   If [], it means that the file did not exist (sts = -1 or -2) or input
%   filename was incorrect (sts = -3).
%
%   f_dlt: If a file is generated in this function, in particular, a .nii
%   from a .nii.gz, it will be also return as f_dlt, i.e., where dlt
%   denotes delete; f_dlt can be safely deleted when needed ensuring
%   origial file is retained. If f_dlt is [], no new file was generated in
%   this function, thus, nothing to cleanup later.
%
%   sts: 1, -1, -2, or -3. If sts=1, all good, that is, file existed and a
%   loadale version has been returned. If sts = -1 or -2, input file or its
%   alternative version (.nii.gz if input was .nii and vice versa) did not
%   exist. Is sts = -3, input file format is not recongnized; teh oyther
%   two outputs are retunred as [].
%
% Hamid Behjat

sts = 1;
if endsWith(f, '.nii')
    if exist(f, 'file')
        f_dlt = [];
    else
        f_gz = [f, '.gz'];
        if exist(f_gz, 'file')
            gunzip(f_gz);
            f_dlt = f;
        else
            sts = -1; % missing input .nii file, or its .nii.gz
        end
    end
elseif endsWith(f, '.nii.gz')
    f_gz = f;
    f = strrep(f_gz, '.nii.gz', '.nii');
    if exist(f, 'file')
        f_dlt = [];
    else
        if exist(f_gz, 'file')
            gunzip(f_gz);
            f_dlt = f;
        else
            sts = -2; % missing input .nii.gz file, or its .nii
        end
    end
else
    sts = -3; % unknown file
end
if sts~=1
    f = [];
    f_dlt = [];
end
end