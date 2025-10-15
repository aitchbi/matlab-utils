function N = hb_annot_get_parcsize(f_p)
[vtx_p, s_p, tbl_p] = read_annotation(f_p);
assert(isequal(length(vtx_p),length(s_p)));
clbls = tbl_p.table(:,5); % color-coded lbls
assert(isequal(sort(clbls), sort(unique(s_p))));
assert(length(clbls)==tbl_p.numEntries);
I_bg = 1;
lbl_bg = 65793; % background; RGB:[1 1 1]
assert(tbl_p.table(I_bg,5)==lbl_bg, 'Background label not found');
assert(contains(tbl_p.struct_names{I_bg}, 'Background'));
assert(contains(tbl_p.struct_names{I_bg}, 'Medial_Wall'));
Nparc = tbl_p.numEntries-1;
clbls = clbls(2:end);
N = zeros(Nparc,1);
for iP=1:Nparc
    N(iP) = nnz(s_p==clbls(iP));
end
end