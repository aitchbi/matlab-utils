function hb_progress(n,N,varargin)
% HB_PROGRESS displays progress.
%
% h behjat 

d = inputParser;
addParameter(d,'init',0);
addParameter(d,'tag','');
addParameter(d,'l',length(num2str(N)));
parse(d,varargin{:});
opts = d.Results;

if n==1 || opts.init
    if isempty(opts.tag)
        fprintf(['\n','progress: ']);
    else
        fprintf(['\n',opts.tag]);
    end
else
    fprintf(repmat('\b',1,2*opts.l+1),n);
end
l = num2str(opts.l);
eval(['fprintf(''%-',l,'d/%-',l,'d'',n,N)']);