function [lbls, s_p, clbls, BG] = hb_gii_get_parc_labels(f_parc)
% HB_GII_GET_PARC_LABELS reads input surface parcellation file (eg DK or
% Schaefer) and returns labels and other info. 
%
% 
% h behjat 

[vtx_p, s_p, tbl_p] = read_annotation(f_parc);

assert(isequal(length(vtx_p),length(s_p)));

clbls = tbl_p.table(:,5); % color-coded lbls
assert(isequal(sort(clbls), sort(unique(s_p))));
assert(length(clbls)==tbl_p.numEntries);

I_bg = 1;

lbl_bg = 65793; % background; RGB:[1 1 1]

assert(tbl_p.table(I_bg,5)==lbl_bg, 'Background label not found');

if contains(f_parc, 'Yeo2011_7Networks_N1000')

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
