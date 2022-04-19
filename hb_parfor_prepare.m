function p = hb_parfor_prepare(L,varargin)
%HB_PARFOR_PREPARE initiates a parallel pool. The number of workers will be
%the maximum number of workers available on the local profile, unless L
%(length of for loop) is smaller than that. In the latter case, the size of
%the pool will be equal to L. If second input is provided, i.e., desired
%size of pool, the size of the initiated pool will be equal or smaller to
%that value.
%
% Inputs: 
%   L: length of for loop to run in parallel.
%   varargin{1}: desired parpool size.
%
% Hamid Behjat

if nargin<2
    Np = [];
else
    Np = varargin{1};
end

d = parcluster('local'); 
Nmax = d.NumWorkers; % max possible # of workers

if isempty(Np)
    Np = Nmax;
else
    if Np>Nmax
        s = ['[HB] desired parpool size ',...
            'larger than max possible pool size;',...
            ' a smaller pool is initiated.'];
        warning(s);
        Np = Nmax;
    end
end

if L<Np
    Np=L;
end

p = gcp('nocreate');
if isempty(p)
    p = parpool(Np);
elseif p.NumWorkers~=Np
    delete(p);
    p = parpool(Np);
end
end


