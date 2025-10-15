function opts = hb_get_atlasinfo(WhichAtlas, d_fsaverage_hb, d_fsmatlab)

if not(exist('read_annotation.m','file'))
    if exist('d_fsmatlab', 'var')
        addpath(d_fsmatlab); % read_annotation.m
    else
        error('read_annotation.m or freesurfer matlab folder not found.');
    end
end

Yeo7 = {
    'Vis'
    'SomMot'
    'DorsAttn'
    'SalVentAttn'
    'Limbic'
    'Cont'
    'Default'};

Yeo7_v2 = {
    'VN'
    'SMN'
    'DAN'
    'VAN'
    'LN'
    'FPN'
    'DMN'};

SchaeferYeo7 = {
    'Schaefer100Yeo7'
    'Schaefer200Yeo7'
    'Schaefer300Yeo7'
    'Schaefer400Yeo7'
    'Schaefer500Yeo7'
    'Schaefer600Yeo7'
    'Schaefer700Yeo7'
    'Schaefer800Yeo7'
    'Schaefer900Yeo7'
    'Schaefer1000Yeo7'
    };

SchaeferYeo17 = {
    'Schaefer100Yeo17'
    'Schaefer200Yeo17'
    'Schaefer300Yeo17'
    'Schaefer400Yeo17'
    'Schaefer500Yeo17'
    'Schaefer600Yeo17'
    'Schaefer700Yeo17'
    'Schaefer800Yeo17'
    'Schaefer900Yeo17'
    'Schaefer1000Yeo17'
    };

SchaeferYeo = [SchaeferYeo7; SchaeferYeo17];

getSchYeo7 = @(x) sprintf('Schaefer2018_%dParcels_7Networks_order',x);

getSchYeo17 = @(x) sprintf('Schaefer2018_%dParcels_17Networks_order',x);

d_lbls = fullfile(d_fsaverage_hb,'label');
switch WhichAtlas
    case 'Glasser360'
        d = 'HCPMMP1';
        opts.Nroi_hemi = 180;
        
    case 'Yeo17'
        d = 'Yeo2011_17Networks_N1000';
        opts.Nroi_hemi = 17;
        
    case 'Yeo7'
        d = 'Yeo2011_7Networks_N1000';
        opts.Nroi_hemi = 7;
        
    case 'Yeo17_volumetric'
        opts.Nroi_hemi = 'n/a'; % not fixed; depends on MinRoiSize
        
    case SchaeferYeo
        d_str = 9;
        d_end = strfind(WhichAtlas,'Yeo')-1;
        opts.Nroi_hemi = str2double(WhichAtlas(d_str:d_end))/2;
        
        d = opts.Nroi_hemi*2;
        
        switch WhichAtlas
            case SchaeferYeo7
                d = getSchYeo7(d);

            case SchaeferYeo17
                d = getSchYeo17(d);
        end
        opts.Nroi = opts.Nroi_hemi*2;
    
    case {'DesikanKilliany', 'DK'}
        d = 'aparc';
        opts.Nroi_hemi = 35;
        opts.aparc_label_order.lh  = 1:opts.Nroi_hemi;
        opts.aparc_label_order.rh  = 1:opts.Nroi_hemi;
        opts.aparc_label_number.lh = 1000 + opts.aparc_label_order.lh;
        opts.aparc_label_number.rh = 2000 + opts.aparc_label_order.rh;

    case 'Braak3'
        opts.Nroi_hemi = 3;
        d = [];

    case 'Braak5'
        opts.Nroi_hemi = 5;
        d = [];
end
opts.ParcName = d;

switch WhichAtlas
    case {'Braak3', 'Braak5'}
        % no annotation file

    otherwise
        f_lh = fullfile(d_lbls,sprintf('lh.%s.annot',opts.ParcName));
        f_rh = fullfile(d_lbls,sprintf('rh.%s.annot',opts.ParcName));

        [~,~,coltb_lh] = read_annotation(f_lh);
        [~,~,coltb_rh] = read_annotation(f_rh);
end

switch WhichAtlas
    case {'DesikanKilliany', 'DK'}
        
        assert(coltb_lh.numEntries==(1 + opts.Nroi_hemi), 'fishy'); % one is the "unknown" label (medial wall)
        assert(coltb_rh.numEntries==(1 + opts.Nroi_hemi), 'fishy'); % "

        rows_lh = 1 + opts.aparc_label_order.lh;
        rows_rh = 1 + opts.aparc_label_order.rh;

        opts.n_seed.lh = coltb_lh.struct_names(rows_lh); % region names
        opts.n_seed.rh = coltb_rh.struct_names(rows_rh);

        opts.aparc_label_color.lh = coltb_lh.table(rows_lh, 5)'; % region color-codes
        opts.aparc_label_color.rh = coltb_rh.table(rows_rh, 5)';

    case 'Braak3'
        opts.n_seed.lh = {
            'Braak-I-II'
            'Braak-III-IV'
            'Braak-V-VI'
            };
        opts.n_seed.rh = opts.n_seed.lh;

    case 'Braak5'
        opts.n_seed.lh = {
            'Braak-I-II'
            'Braak-III'
            'Braak-IV'
            'Braak-V'
            'Braak-VI'
            };
        opts.n_seed.rh = opts.n_seed.lh;

    otherwise
        opts.n_seed.lh = coltb_lh.struct_names(2:end);
        opts.n_seed.rh = coltb_rh.struct_names(2:end);
end

for iSeed=1:opts.Nroi_hemi
    
    d = opts.n_seed.lh{iSeed};
    switch WhichAtlas
        case 'Glasser360'
            d = strrep(d,'L_','');
            d = strrep(d,'_ROI','');
            
        case 'Yeo17'
            d = strrep(d,'17Networks_','');
            
        case 'Yeo7'
            d = strrep(d,'7Networks_','');
            
        case SchaeferYeo7
            d = strrep(d,'7Networks_LH_','');
            d = strrep(d,'_',' ');
            [inetw, netw] = getnetw(d, Yeo7);
            opts.n_netw.lh{iSeed} = netw;
            opts.n_netwnum.lh(iSeed) = inetw;
            
        case SchaeferYeo17
            d = strrep(d,'17Networks_LH_','');
            d = strrep(d,'_',' ');
        
        case {'DesikanKilliany', 'DK', 'Braak3', 'Braak5'}
            % no change
    end
    opts.n_seed.lh{iSeed} = d;
    
    d = opts.n_seed.rh{iSeed};
    switch WhichAtlas
        case 'Glasser360'
            d = strrep(d,'R_','');
            d = strrep(d,'_ROI','');
            
        case 'Yeo17'
            d = strrep(d,'17Networks_','');
            
        case 'Yeo7'
            d = strrep(d,'7Networks_','');
            
        case SchaeferYeo7
            d = strrep(d,'7Networks_RH_','');
            d = strrep(d,'_',' ');
            [inetw, netw] = getnetw(d, Yeo7);
            opts.n_netw.rh{iSeed} = netw;
            opts.n_netwnum.rh(iSeed) = inetw;

        case SchaeferYeo17
            d = strrep(d,'17Networks_RH_','');
            d = strrep(d,'_',' ');
        
        case {'DesikanKilliany', 'DK', 'Braak3', 'Braak5'}
            % no change
    end
    opts.n_seed.rh{iSeed} = d;
end

switch WhichAtlas
    case SchaeferYeo7
        opts.Yeo7    = Yeo7;
        opts.Yeo7_v2 = Yeo7_v2;

        opts = get_ntwinfo(opts);

        opts.n_netwnum_cell.lh   = zeros(7, opts.Nroi_hemi);
        opts.n_netwnum_cell.rh   = zeros(7, opts.Nroi_hemi);
        for iNetw=1:7
            opts.n_netwnum_cell.lh(iNetw, :) = opts.n_netwnum.lh==iNetw;
            opts.n_netwnum_cell.rh(iNetw, :) = opts.n_netwnum.rh==iNetw;
        end
    
    otherwise
        % extend if needed
end

%fprintf('\n\n..NOTE: Assuming 1st half of FC is lh and 2nd half is rh.\n\n');
end

%==========================================================================
function [netwnum, netw] = getnetw(d,yeo7)
netwnum = find(cellfun(@(x) contains(d, x), yeo7));
netw = yeo7{netwnum};
end

%==========================================================================
function opts = get_ntwinfo(opts)

ntwinfo = struct;

ntwinfo.I = [
    opts.n_netwnum.lh(:)
    opts.n_netwnum.rh(:)
    ];

N_ntw = length(unique(ntwinfo.I));

ntwinfo.L_ntw = zeros(N_ntw,1);
ntwinfo.I_ntw_cell = cell(N_ntw,1); % [29.10.2024]
ntwinfo.I_sort = [];
for k=1:N_ntw
    d = find(ntwinfo.I==k);
    ntwinfo.I_ntw_cell{k} = d;
    ntwinfo.L_ntw(k) = length(d);
    ntwinfo.I_sort = [ntwinfo.I_sort d(:)'];
end

ntwinfo.I_ntwsort = ntwinfo.I(ntwinfo.I_sort);

switch N_ntw
    case 7
        ntwinfo.minilabels = cell(1,N_ntw);
        ntwinfo.minilabels{1} = 'Visual';
        ntwinfo.minilabels{2} = 'Motor';
        ntwinfo.minilabels{3} = 'DAN';
        ntwinfo.minilabels{4} = 'VAN';
        ntwinfo.minilabels{5} = 'Limbic';
        ntwinfo.minilabels{6} = 'FPCN';
        ntwinfo.minilabels{7} = 'DMN';
    otherwise
end

ntwinfo.TickNetw = cumsum(ntwinfo.L_ntw)-floor(ntwinfo.L_ntw/2);


ntwinfo.I_ntwsort_cell = cell(1,7);
for k=1:7
    d = find(ntwinfo.I_ntwsort==k,1);
    ntwinfo.I_ntwsort_cell{k} = d:d+ntwinfo.L_ntw(k)-1; % [29.20.2024] changed from "ntwinfo.I_ntw" to "ntwinfo.I_ntwsort_cell" as otheise its prone to be misinterpretted as unsorted indices. Instead, a "ntwinfo.I_ntw_cell" field now also exists.
end

opts.netwinfo = ntwinfo;
end
