function [S_dcomp,S_dcomp_full] = hb_sosks_decompose(S,Nsbs,L,c_sosks,lmax,varargin)
%HB_SOSKS_DECOMPOSE Decomposing a given set of graph signals S using a
% system of spectral kernels defined on the graph Laplacian spectrum using
% Chebyshev polynomila approximations of the kernels (given by their
% Chebyshev polynomial coefficients c_sosks).
%
% Requires: 
% sgwt_cheby_op.m etc. 
%
% Hamid Behjat 

p = inputParser;
addParameter(p,'sbs_todo',[]);
addParameter(p,'alsoSaveFullDecomp',false); % also save full decomposed signals?
addParameter(p,'f_dcomp',[]);      % check whether file partially computed & saved before; see NOTE 1.
addParameter(p,'f_dcomp_full',[]); % check whether file partially computed & saved before; ; see NOTE 1.
addParameter(p,'ss',true); % verbose?
parse(p,varargin{:});
opts = p.Results;

M=50;

if isempty(opts.sbs_todo)
    subbands = 1:Nsbs;
else
    subbands = find(opts.sbs_todo);
end

Nss = size(S,2);

if opts.alsoSaveFullDecomp
    Ng = length(L);
end

if ~isempty(opts.f_dcomp)
    [p,n,e] = fileparts(opts.f_dcomp);
    f_dcomp_partial = fullfile(p,[n,'_partialRun',e]);
else
    f_dcomp_partial = [];
end

if ~isempty(opts.f_dcomp_full)
    [p,n,e] = fileparts(opts.f_dcomp_full);
    f_dcomp_full_partial = fullfile(p,[n,'_partialRun',e]);
else
    f_dcomp_full_partial = [];
end

if exist('f_dcomp_partial','file')
    d = load(f_dcomp_partial);
    S_dcomp = d.S_dcomp;
    d = size(S_dcomp,2);
    S_dcomp = cat(2,S_dcomp,zeros(Nsbs,Nss-d));
    iF = d+1; %#ok<NASGU>
else
    S_dcomp = zeros(Nsbs,Nss);
    iF = 1; %#ok<NASGU>
end

if exist('f_dcomp_full_partial','file')
    d = load(f_dcomp_full_partial);
    S_dcomp_full = d.S_dcomp_full;
    d = size(S_dcomp_full,2);
    S_dcomp_full = cat(2,S_dcomp_full,zeros(Nsbs,Nss-d,Ng));
    iF = d+1; % see NOTE 2.
else
    S_dcomp_full = zeros(Nsbs,Nss,Ng);
    iF = 1;
end

skipPartialSave = 0;

for k=iF:Nss
    d = sgwt_cheby_op(S(:,k),L,c_sosks,[0,lmax]);
    for iB = subbands
        S_dcomp(iB,k) = norm(real(d{iB}))^2;
        
        if opts.alsoSaveFullDecomp
            S_dcomp_full(iB,k,:) = d{iB}(:);
        end
    end
    if ~rem(k,M)
        if ~isempty(f_dcomp_partial)
            d = S_dcomp;
            S_dcomp = S_dcomp(:,1:k);
            save(f_dcomp_partial,'S_dcomp');
            S_dcomp = d;
        end
        if ~isempty(f_dcomp_full_partial) && ~skipPartialSave
            d = S_dcomp_full;
            S_dcomp_full = S_dcomp_full(:,1:k,:);
            d1 = whos('S_dcomp_full');
            if d1.bytes<1.9e9
                save(f_dcomp_full_partial,'S_dcomp_full');
            else
                skipPartialSave = 1;
                warning('HB: large file! partial run save skipped');
                % see NOTE 3.
            end
            S_dcomp_full = d;
        end
    end
    if opts.ss
        prgs(k,Nss,0,'# of frames decomposed:');
    end
end

if ~opts.alsoSaveFullDecomp
    S_dcomp_full = [];
end

end

function prgs(n,N,init,tag)
l = numel(num2str(N));
if n==1 || init
    fprintf(['\n',tag]);
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end

% NOTE 1:
% when the run is long, for instrance when computing 'allFrames', it's good
% to so save partial runs, 50 frames at a time, in case you need to break
% or it times out (lunarc!). Next time you run, loads previous saves and
% continues from there; this is especially good on lunarc, or even on golgi
% etc. if you happen to need to break in the middle of a run. 
%
% NOTE 2:
% Note that this iF is overwriting the one above. But that's ok, since most
% likely it should be the same as the one above. However, it may be smaller
% than the one above, in the case that the run was interuppted just at the
% time when S_dcomp_full_partialRun was being save. In this case, we would
% want to rerun from this smaller iF.
%
% NOTE 3:
% Partial run saves were incorporated since it typically happened that
% parallel runs on golgi were crashed, and thus, it wa good to do partial
% run saves to save computational time in consequent runs. But doing
% partial run save for f_dcomp_full doesn't make much sense of the file is
% super large. f_dcomp_full is intended for computing either:
%
% i) only few subbands, like
% lowpass/highpass to compute strutural decoupling index. For instance, for
% a 100K node graph, for ~ 300 signals, it will be around 0.5 GB. If you
% want to save larger results (probably not), adjust save commands to use
% v7.3.
%
% ii) [aug 2021] computing all subbands, for running k-means clustering on the energy
% profile of all voxels in a white matter hemipher; in this case, the file
% will be uper large, and it just does not make sense to save it, not only
% for space reasons but also for the time it takes. In this case, The file
% that is generated should just be output by the function, be processed,
% and then discarded. 

