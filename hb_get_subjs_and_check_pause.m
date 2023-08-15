function [I, PauseLength] = hb_get_subjs_and_check_pause(n_file, RunPause)

% Inputs:
%   n_file: file name of a script that is to be run, without .m extension. 
% 
% Example usuage: 
% 
% [I, PauseLength] = hb_get_subjs_and_pause(mfilename);
% 
% *** Expected format of n_file:
%
% option 1: '*_SubjsFrom<X>To<Y>'
%
% option 2: '*_SubjsFrom<X>To<Y>_Pause<T1><T2>'
%
%  <X>: an integer; first subject number 
%  <Y>: an integer; last subject number
% <T1>: an integer; pause length (minutes or hours)
% <T2>: 'Minutes' or 'Hours'
% 
% Outputs:
%   I: subject numbers
%   PauseLength: in seconds; to feed into function pause.m
%
% Hamid Behjat

if not(exist('RunPause','var')) || isempty(RunPause)
    RunPause = false;
end

d1 = strfind(n_file, 'From') + length('From');
d2 = strfind(n_file, 'To');
d3 = d2 + length('To');
iSta = str2double(n_file(d1:d2-1));

if contains(n_file, 'Pause')
    iPause = strfind(n_file, 'Pause') + length('Pause');
    if contains(n_file, 'Hours')
        iHours = strfind(n_file, 'Hours');
        PauseLength = str2double(n_file(iPause:iHours-1));
        PauseLength = PauseLength*60*60; % hours to seconds
    else
        assert(contains(n_file, 'Minutes'),...
            'Pause should be give in either Hours or Minutes.');
        iMinutes = strfind(n_file, 'Minutes');
        PauseLength = str2double(n_file(iPause:iMinutes-1));
        PauseLength = PauseLength*60; % minutes to seconds
    end
    iEnd = str2double(n_file(d3:strfind(n_file, '_Pause')-1));
else
    PauseLength = 0;
    iEnd = str2double(n_file(d3:end));
end

I = iSta:iEnd;

if RunPause
    if PauseLength~=0
        fprintf('\n.Script paused for %d seconds at %s. \n', PauseLength, datetime);
        pause(PauseLength);
    end
end
end