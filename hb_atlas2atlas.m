function [y, yy, xlhrh] = hb_atlas2atlas(x, atlas_x, atlas_y, d_fsavg, varargin)
%
% y: parcel values for atlas_y
% yy: surface vertex values associated to each parcel in atlas_y
% xlhrh: surface vertex map
%
% h behjat

d = inputParser;
addParameter(d,'Hemisphere', []);
parse(d,varargin{:});
opts = d.Results;

[x, Nx] = xchk(x);

d = hb_get_atlasinfo_mini(atlas_x);
Nr = d.Nr;
Nr_hemi = d.Nr_hemi;

if Nx==Nr

    hemis = {'lh', 'rh'};
    ylhrh = struct;
    xlhrh = struct;
    for iH=1:2
        hemi = hemis{iH};
        switch hemi
            case 'lh'
                xhemi = x(1:Nr_hemi);
            case 'rh'
                xhemi = x(Nr_hemi+1:end);
        end
        [ylhrh.(hemi), yylhrh.(hemi), xlhrh.(hemi)] = runpvt(xhemi, atlas_x, atlas_y, hemi, d_fsavg);
    end
    y = ylhrh2y(ylhrh);
    yy = yylhrh2y(yylhrh);

elseif Nx==Nr_hemi

    hemi = opts.Hemisphere;
    assert(~isempty(hemi));
    [y, yy] = runpvt(x, atlas_x, atlas_y, hemi, d_fsavg);

else
    error('fishy');
end
end

%==========================================================================
function [y_mean, yy, xx] = runpvt(x, atlas_x, atlas_y, hemi, d_fsavg)
[~, lbls_x, N_hemi_x] = hb_get_surfatlas(atlas_x, hemi, d_fsavg, 'white');
[~, lbls_y, N_hemi_y] = hb_get_surfatlas(atlas_y, hemi, d_fsavg, 'white');
assert(isequal(length(lbls_x),length(lbls_y)));
ulbls_x = unique(lbls_x);
ulbls_y = unique(lbls_y);
assert(ulbls_x(1)==0);
assert(ulbls_y(1)==0);
ulbls_x = ulbls_x(2:end);
ulbls_y = ulbls_y(2:end);
L_x = length(ulbls_x);
L_y = length(ulbls_y);
assert(isequal(L_x, N_hemi_x));
assert(isequal(L_y, N_hemi_y));
N_surf = length(lbls_x);
xx = zeros(N_surf,1);
for k=1:L_x
    lbl = ulbls_x(k);
    xx(lbls_x==lbl) = x(lbl);
end
y_mean = zeros(L_y,1);
yy = cell(L_y,1);
for k=1:L_y
    lbl = ulbls_y(k);
    I = lbls_y==lbl;
    y_mean(lbl) = mean(xx(I));
    yy{lbl} = xx(I);
end
end

%==========================================================================
function y = ylhrh2y(x)
Nlh = length(x.lh);
Nrh = length(x.rh);
N = Nlh + Nrh;
y = zeros(N,1);
y(1:Nlh)     = x.lh;
y(Nlh+1:end) = x.rh;
end

%==========================================================================
function yy = yylhrh2y(x)
Nlh = length(x.lh);
Nrh = length(x.rh);
N = Nlh + Nrh;
yy = cell(N,1);
yy(1:Nlh)     = x.lh;
yy(Nlh+1:end) = x.rh;
end

%==========================================================================
function [x, Nx] = xchk(x)
assert(isvector(x));
x = x(:);
Nx = length(x);
end
