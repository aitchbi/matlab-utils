function [surf,labels,Nr_hemi,lbls_unique] = hb_get_surfatlas(atlas,hemi,FSavgDir,WhichSurf,varargin)
%
%
%
% 
% h behjat 

d = inputParser;
addParameter(d,'SkipAssertionChecks', 0);
parse(d,varargin{:});
opts = d.Results;

switch atlas
    case 'Braak5'
        n_atlas = 'aparc';
        Nr_hemi = 5;
        lbls = 1:5;
        % 1: I-II
        % 2: III
        % 3: IV
        % 4: V
        % 5: VI
    
    case 'Braak3'
        n_atlas = 'aparc';
        Nr_hemi = 3;
        lbls = 1:3;
        % 1: I-II
        % 2: III-IV
        % 3: V-VI

    case 'EBM5'
        n_atlas = 'aparc';
        Nr_hemi = 5;
    
    otherwise
        param = hb_get_atlasinfo_mini(atlas,hemi);
        n_atlas = param.ParcName;
        Nr_hemi = param.Nr_hemi;
        lbls    = param.lbls;
end

f_surf = fullfile(FSavgDir,'surf',[hemi,'.',WhichSurf]);
f_anot = fullfile(FSavgDir,'label',[hemi,'.',n_atlas,'.annot']);

% surface
surf = struct;
[surf.vertices, surf.faces, ~] = read_surf(f_surf);
surf.faces = surf.faces + 1;

% atlas

IsNonZeroLabelForUnkownOrUnassignedVerticesUsed = false;

switch atlas
    case {'Braak5', 'Braak3', 'EBM5'}


        LabelForUnkownOrUnassignedVertices = 6;
        % If set to 0, theses regions will have same color code as medial wall, but
        % note that the parcel bounaries will get somewhat warped around these
        % regions.
        %
        % If you assign to a number label, that will prevent messing up parcel
        % boundaries.

        [BRK, EBM] = hb_braak_info([], 'JustgetLabels', 1);
        BRKL = BRK.Labels;
        EBML = EBM.Labels;

        [~, lbls_aparc, tbl_aparc] = read_annotation(f_anot);

        ainfo_DK = hb_get_atlasinfo('DK', FSavgDir);
        
        label_orders  = ainfo_DK.aparc_label_order.(hemi);
        label_numbers = ainfo_DK.aparc_label_number.(hemi);
        label_colors  = ainfo_DK.aparc_label_color.(hemi);

        verifylbls(lbls_aparc, tbl_aparc, ainfo_DK, hemi);
        
        switch hemi
            case 'lh'
                hemi_base = 1000;
            case 'rh'
                hemi_base = 2000;
        end
        
        labels = zeros(size(lbls_aparc));

        for iR = [0 label_orders]

            if iR==0
                lbliR_color  = tbl_aparc.table(1,5);
                lbliR_number = hemi_base;
            else
                lbliR_color  = label_colors(iR);
                lbliR_number = label_numbers(iR);
            end

            switch atlas
                case 'Braak5'
                    switch lbliR_number
                        case {1000 2000} % unkown vertices
                            lbl = LabelForUnkownOrUnassignedVertices;
                        case num2cell(BRKL.I_II)
                            lbl = 1; % Braak I-II
                        case num2cell(BRKL.III)
                            lbl = 2; % Braak III
                        case num2cell(BRKL.IV)
                            lbl = 3; % Braak IV
                        case num2cell(BRKL.V)
                            lbl = 4; % Braak V
                        case num2cell(BRKL.VI)
                            lbl = 5; % Braak VI
                        case {1004, 2004} % Corpus Callosum
                            d = lbls_aparc==lbliR_color;
                            assert(nnz(d)==0, ...
                                'fishy: corpus callosum not expected');
                            continue;
                        otherwise
                            d1 = ismember(lbliR_number, BRKL.between_IV_V);
                            d2 = ismember(lbliR_number, BRKL.between_V_VI);
                            assert(or(d1,d2), 'fishy: unknown aparc label');
                            %lbl = 0; % same color as unknown 
                            lbl = LabelForUnkownOrUnassignedVertices; 
                    end

                case 'Braak3'
                    switch lbliR_number
                        case {1000 2000}
                            lbl = LabelForUnkownOrUnassignedVertices;
                        case num2cell(BRKL.I_II)
                            lbl = 1; % Braak I-II
                        case num2cell(BRKL.composite_III_IV)
                            lbl = 2; % Braak III-IV
                        case num2cell(BRKL.composite_V_VI)
                            lbl = 3; % Braak V-VI
                    end
            
                case 'EBM5'
                    switch lbliR_number
                        case {1000 2000}
                            lbl = LabelForUnkownOrUnassignedVertices;
                        case num2cell(EBML.I)
                            lbl = 1; % EBM_I
                        case num2cell(EBML.II)
                            lbl = 2; % EBM_II
                        case num2cell(EBML.III)
                            lbl = 3; % EBM_III
                        case num2cell(EBML.IV)
                            lbl = 4; % EBM_IV
                        case num2cell(EBML.V)
                            lbl = 5; % EBM_V
                    end
            end
            
            if 0 % debug
                if lbl==LabelForUnkownOrUnassignedVertices
                    fprintf('\n. %s -- label: %d -- number in annotation: %d', ...
                        f_anot, ...
                        lbliR_number, ...
                        nnz(lbls_aparc==lbliR_color));
                end
            end
            
            if LabelForUnkownOrUnassignedVertices~=0
                if nnz(lbls_aparc==lbliR_color)>0
                    IsNonZeroLabelForUnkownOrUnassignedVerticesUsed = true;
                end
            end
            
            labels(lbls_aparc==lbliR_color) = lbl;
        end

        if IsNonZeroLabelForUnkownOrUnassignedVerticesUsed
            ExtraLabel = 1;
        else
            ExtraLabel = 0;
        end

        N_total = Nr_hemi+1;

        switch atlas
            case {'Braak5', 'Braak3', 'EBM5'}
                N_total = N_total+ExtraLabel;
        end

        assert(length(unique(labels))==N_total);

    otherwise

        [~, labels, tbl_aparc] = read_annotation(f_anot);
        for iR = [0 lbls] %0:Nr_hemi
            labels(labels==tbl_aparc.table(iR+1,5)) = iR;
        end
        
        N_total = Nr_hemi+1;
end

lbls_unique = unique(labels);

if opts.SkipAssertionChecks
    return;
else
    d = sort(lbls_unique);
    if ~isequal(length(d), N_total)
        switch atlas
            case {'DesikanKilliany', 'DK'}
                % label 4 (corpuscallosum) not in the cerebral cortex xh.aparc
                if IsNonZeroLabelForUnkownOrUnassignedVerticesUsed
                    dd = [d(:); LabelForUnkownOrUnassignedVertices];
                    dd = sort(dd);
                    ddd = sort([0:3, 5:35, LabelForUnkownOrUnassignedVertices]);
                    assert(isequal(dd(:), ddd(:)));
                else
                    assert(isequal(d(:)', [0:3, 5:35]));
                end
                fprintf('\n.NOTE: corpuscallosum (label 1004/2004) not in xh.aparc');
            otherwise
                error('extend');
        end
    end
    if d(1)~=0
        switch atlas
            case 'Glasser360'
                % In Glasser's atlas, label for medial wall does not match that
                % given in the colortable, so we fix.
                assert(d(end-1)==Nr_hemi);
                assert(~ismember(d(end),1:Nr_hemi));
                labels(labels==d(end)) = 0;
            otherwise
                error('extend');
        end
    end
end
end

%==========================================================================
function verifylbls(lbls_aparc, tbl_aparc, ainfo_DK, hemi)
label_colors = ainfo_DK.aparc_label_color.(hemi);
ucolors = unique(lbls_aparc);
assert(isequal(tbl_aparc.struct_names{1}, 'unknown')); % label for "unkown" region
lbl_unkown = tbl_aparc.table(1,5);
label_ok_to_be_missing = [0 lbl_unkown]; % ok labels to be missing in annotation file; 0 is unknown
label_missing = setdiff(ucolors,label_colors);
chk1 = isempty(label_missing);
chk2 = all(ismember(label_missing, label_ok_to_be_missing));
if or(chk1, chk2)
    return;
else
    chk1 = length(ucolors);
    chk2 = ainfo_DK.Nroi_hemi+1;
    assert(chk2>chk1, ...
        'fishy: more labels in annotation file than labels in label table');
    for k=1:length(label_missing)
        switch label_missing(k)
            case lbl_unkown
                % OK; ok not to have any "unknown" labels in the annotation
            case 0
                % OK; background i.e. medial wall
            otherwise
                d = label_missing(k)==tbl_aparc.table(:,5);
                error('Fishy: label in annotation missing in label table: %s', ...
                    tbl_aparc.struct_names(d));
        end
    end
end
end
