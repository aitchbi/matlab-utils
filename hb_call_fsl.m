function [sts,log] = hb_call_fsl(cmd,fsldir)
%HB_CALL_FSL initiates FSL, then calls the requested command.
%
% inputs:
%   cmd: FSL command as you would input in Terminal. 
%
%   fsldir: if FSL is already installed, this input is not needed,
%   otherwise, specify.
%
% examples:
% hb_call_fsl('fast -v /absolute/file/path/t1w.nii');
% hb_call_fsl('fast -v /a/b/c/t1w.nii','/usr/local/fsl');
%
% h behjat

if ~exist('fsldir','var')
    fsldir = [];
end

d1 = 'PATH=${FSLDIR}/bin:${PATH}';
d2 = 'export FSLDIR PATH';
d3 = '. ${FSLDIR}/etc/fslconf/fsl.sh';

[sts,log] = system('$FSLDIR');
if sts==126 && contains(log,': is a directory')
    cmd = sprintf('%s; %s; %s; %s',d1,d2,d3,cmd);
else
    assert(~isempty(fsldir),'FSLDIR not set; specify as input.')
    d0 = sprintf('FSLDIR=%s',fsldir);
    cmd = sprintf('%s; %s; %s; %s; %s',d0,d1,d2,d3,cmd);
end

[sts,log] = system(cmd);
end

