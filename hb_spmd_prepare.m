function [y,p] = hb_spmd_prepare(x,splitwhich,varargin)
% HB_SPMD_PREPARE prepares inputs and parallel pool for an spmd run. 
%
% INPUTS:
%   x: cell array of all inputs necessary for a function to be run in spmd.
% 
%   splitwhich: a logical array of size x, specifying which elements in x 
%   are to be split across workers. OK to use 0/1 instead of logicals. 
%
%   Name,Value pair arguments:
%   'Np': desired size of parpool.
% 
% OUTPUTS:
%   y: cell array, equal to x in elements that needed not to be split, but 
%   those that had to be split for smpd, have been split.
%
%   p : handle to the parpool which spmd will use. 
%
% Examples:
%   [1]
%   x{1} = var1;
%   x{2} = var2;
%   x{3} = var3;
%   splitwhich = false(size(x)); 
%   splitwhich(1) = true;
%   [y,p] = hb_spmd_prepare(x,splitwhich);
%   spmd
%   function(y{1}{labindex},y{2},y{3});
%   end
%   % withount spmd the function will be called as:   
%   function(x{1},x{2},x{3});
%
%   [2]
%   delete(gcp('nocreate'));
%   [y,p] = hb_spmd_prepare(x,splitwhich,'Np',3);
%   spmd
%   function(y{1}{labindex},y{2},y{3});
%   end
%
% Hamid Behjat

d = inputParser;
addParameter(d,'Np',[]);
parse(d,varargin{:});
Np = d.Results.Np;

if any(~logical(splitwhich))
    assert(all(xor(splitwhich==1,splitwhich==0)));
    splitwhich = logical(splitwhich);
end

d = parcluster('local');
Nmax = d.NumWorkers;

if isempty(Np)
    Np = Nmax;
else
    if Np>Nmax
        Np = Nmax;
        s = ['\n..desired parpool size ',...
            'larger than max possible pool size;',...
            ' a smaller parpool is created.'];
        fprintf(s);
    end
end

inputs_split = x(splitwhich);
L = length(inputs_split{1});
if L<Np
    Np = L;
end

d = cell(1,nnz(splitwhich));
for iP = 1:length(inputs_split)
    d{iP} = hb_splitarray(inputs_split{iP},Np);
end
y = x;
y(splitwhich) = d;

p = gcp('nocreate');
if isempty(p)
    p = parpool(Np);
else
    if Np~=p.NumWorkers
        s = ['\n..open parpool size not suitable; ',...
            'a new parpool is created.'];
        fprintf(s);
        delete(p);
        p = parpool(Np);
    end
end
end

