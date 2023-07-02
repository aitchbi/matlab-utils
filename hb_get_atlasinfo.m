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
end

f_lh = fullfile(d_lbls,sprintf('lh.%s.annot',d));
f_rh = fullfile(d_lbls,sprintf('rh.%s.annot',d));

[~,~,coltb_lh] = read_annotation(f_lh);
[~,~,coltb_rh] = read_annotation(f_rh);

opts.n_seed.lh = coltb_lh.struct_names(2:end);
opts.n_seed.rh = coltb_rh.struct_names(2:end);

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
            
        case SchaeferYeo17
            d = strrep(d,'17Networks_LH_','');
            d = strrep(d,'_',' ');
    end
    opts.n_seed.lh{iSeed} = d;
    [inetw, netw] = getnetw(d, Yeo7);
    opts.n_netw.lh{iSeed} = netw;
    opts.n_netwnum.lh(iSeed) = inetw;
    
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
            
        case SchaeferYeo17
            d = strrep(d,'17Networks_RH_','');
            d = strrep(d,'_',' ');
    end
    opts.n_seed.rh{iSeed} = d;
    [inetw, netw] = getnetw(d, Yeo7);
    opts.n_netw.rh{iSeed} = netw;
    opts.n_netwnum.rh(iSeed) = inetw;
end
opts.Yeo7    = Yeo7;
opts.Yeo7_v2 = Yeo7_v2;

%fprintf('\n\n..NOTE: Assuming 1st half of FC is lh and 2nd half is rh.\n\n');

opts = get_ntwinfo(opts);

opts.n_netwnum_cell.lh   = zeros(7, opts.Nroi_hemi);
opts.n_netwnum_cell.rh   = zeros(7, opts.Nroi_hemi);
for iNetw=1:7
    opts.n_netwnum_cell.lh(iNetw, :) = opts.n_netwnum.lh==iNetw;
    opts.n_netwnum_cell.rh(iNetw, :) = opts.n_netwnum.rh==iNetw;
end
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

ntwinfo.I_sort = [];
for k=1:N_ntw
    d = find(ntwinfo.I==k);
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

opts.netwinfo = ntwinfo;
end
