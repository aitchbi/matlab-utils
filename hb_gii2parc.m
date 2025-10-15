function [p, lbls, P, BG] = hb_gii2parc(f_i, f_p, varargin)
% HB_GII2PARC parcellates surface map using a surface parcellation/atlas.
% The surface map and parcellation should belong to the same left/right
% hemisphere surface.
%
% Inputs: 
%   f_i: gifti file, containing a set of M left/right surface map(s), or an
%   NxM matrix containing a surface map of length N (N: number of vertices
%   in f_p) in each column.
%
%   f_p: annotation file, specifying parcellation of left/righ cortex.  
%
% Outputs: 
%   p: PxM matrix, each column being a parcelated map.   
%
%   lbls: structure with info about the parcels, their index, name, etc. 
%
%   P: PxM cell, each cell being the map values in a parcel. That is, mean
%   of P{i,j} is equal to p(i,j).
%
%   BG: structure with info about background/medial-wall in the
%   parcellation file, its name, number of vertices etc.
% 
% Examples:
% [p, lbls, BG] = hb_gii2parc(lh.maps.gii, lh.schaefer200yeo7.annot);
%
% Dependencies:
% .FreeSurfer matlab
% .github/aitchbi/matlab-utils
%
% HB

p = inputParser;
addParameter(p,'Something', []);
parse(p,varargin{:});
opts = p.Results; %#ok<NASGU> 

%-Load files.
[maps, s_p, clbls, lbls, BG] = loadfiles(f_i, f_p);

Nvtx       = size(maps, 1);
Nmaps      = size(maps, 2);
Nparc      = length(clbls);
p          = zeros(Nparc, Nmaps);
P          = cell(Nparc, Nmaps); % map values in parcels
BG.mapvals = cell(1, Nmaps);  % map values in bg

for iMap=1:Nmaps
    map = maps(:,iMap);
    showprgs(iMap,Nmaps,'Parcellating surface maps..');
    if iMap==1
        Nnan = chkmap(map, iMap, BG);
    else
        chkmap(map, iMap, [], Nnan);
    end
    BG.mapvals{iMap} = map(s_p==BG.lbl);
    for iP=1:Nparc
        P{iP,iMap} = map(s_p==clbls(iP));
        p(iP,iMap) = mean(P{iP,iMap}, 'omitnan');
    end
    assert(sum(cellfun(@length, P(:,iMap)))==(Nvtx-BG.lbl_count));
end
if Nmaps==1
    t = '';
else
    t = 'each';
end
fprintf('\n.Number of NaNs (likely medial-wall) in %s map: %d', t, Nnan);
fprintf('\n.Number of background+medial-wall vertices in parcellation: %d', BG.lbl_count);
end


%==========================================================================
function showprgs(n,N,tag)
l = numel(num2str(N));
if n==1
    fprintf('\n..%s ',tag);
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end

%==========================================================================
function [maps, s_p, clbls, lbls, BG] = loadfiles(f_i, f_p)
[vtx_p, s_p, tbl_p] = read_annotation(f_p);
if ismatrix(f_i)
    maps = f_i;
else
    d = gifti(f_i);
    maps = d.cdata;
end
assert(isequal(length(vtx_p),length(s_p)));
assert(ismatrix(maps));
assert(isequal(length(vtx_p), size(maps,1))); % NOTE1

clbls = tbl_p.table(:,5); % color-coded lbls
assert(isequal(sort(clbls), sort(unique(s_p))));
assert(length(clbls)==tbl_p.numEntries);

I_bg = 1;
lbl_bg = 65793; % background; RGB:[1 1 1]
assert(tbl_p.table(I_bg,5)==lbl_bg, 'Background label not found');
if contains(f_p, 'Yeo2011_7Networks_N1000')
    assert(contains(tbl_p.struct_names{I_bg}, 'Medial_Wall'));
else
    assert(contains(tbl_p.struct_names{I_bg}, 'Background'));
    assert(contains(tbl_p.struct_names{I_bg}, 'Medial_Wall'));
end
Nparc = tbl_p.numEntries-1;

clbls = clbls(2:end);

lbls       = struct;
lbls.num   = 1:Nparc;    
lbls.names = tbl_p.struct_names(2:end); % assuming I_bg==1

BG = struct; % background/non-parcel
BG.lbl         = lbl_bg;
BG.lbl_indices = s_p==lbl_bg;
BG.lbl_count   = nnz(BG.lbl_indices);
BG.lbl_name    = tbl_p.struct_names{I_bg};
end

%==========================================================================
function Nnan = chkmap(map, iMap, BG, Nnan)
if iMap==1
    Nnan = nnz(isnan(map));
    assert(Nnan<=BG.lbl_count, 'More thans than expected'); % NOTE2
    msg = 'NaNs not in bg/medial-wall of parcellation file';
    assert(all(ismember(isnan(map), BG.lbl_indices)), msg);
else
    assert(nnz(isnan(map))==Nnan, 'NaNs differ across maps');
end
end

% NOTES.
%--------------------------------------------------------------------------
% NOTE1 We check the number of vertices but assume that the "order" of
% vertices in the files match, i.e. rows of vtx_p and maps. This should
% very likely be the case but need to add some sort of check to ensure.
%
% NOTE2 The maps may have NaN values. These are typically at the Medial
% Wall, as imposed by the function used to generate them (e.g.
% hb_vol2surf). The parcellation file also does not have labels in the
% Medial Wall, an dpotentially some other regions; these are denoted as
% "background" vertices. The number of NaNs in the maps (medial wall) is
% not expected to be greater than the number of background vertices (medial
% wall etc.), and this is waht is checked here. If needed,this assertion
% can be relaxed an skipped (i.e. allowing NaNs at vertices that have exact
% labels in the parcellation, i.e., non-background/medial-wall vertices)
% since when computing the average of vertices within a lable NaNs are
% omitted.
