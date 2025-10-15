function param = hb_get_atlasinfo_mini(atlas,hemi)

% hemisphere needed for 'PALS_B12_Brodmann' atals
if ~exist('hemi','var')
    hemi = [];
end

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

param = struct;

switch atlas
    
    case 'Glasser360'
        Nr_hemi = 180;
        ParcName = 'HCPMMP1';
        lbls = 1:180;

    case 'Yeo17'
        Nr_hemi = 17;
        ParcName = 'Yeo2011_17Networks_N1000';
        lbls = 1:17;

    case 'Yeo7'
        Nr_hemi = 7;
        ParcName = 'Yeo2011_7Networks_N1000';
        lbls = 1:7;

    case SchaeferYeo
        d_str = 9;
        d_end = strfind(atlas,'Yeo')-1;
        Nr_hemi = str2double(atlas(d_str:d_end))/2;
        
        d = Nr_hemi*2;
        
        switch atlas
            case SchaeferYeo7
                ParcName = getSchYeo7(d);
                
            case SchaeferYeo17
                ParcName = getSchYeo17(d);
        end
        
        lbls = 1:Nr_hemi;

    case 'PALS_B12_Brodmann'
        ParcName = 'PALS_B12_Brodmann';
        switch hemi
            case 'lh'
                Nr_hemi = 141;
            case 'rh'
                Nr_hemi = 173;
            otherwise
                Nr_hemi = [];
        end
        lbls = []; % to do

    case 'Destrieux148'
        ParcName = 'aparc.a2009s';
        switch hemi
            case 'lh'
                Nr_hemi = 74;
            case 'rh'
                Nr_hemi = 74;
            otherwise
                Nr_hemi = [];
        end
        lbls = []; % to do

    case {'DesikanKilliany', 'DK'}
        ParcName = 'aparc';
        switch hemi
            case 'lh'
                %Nr_hemi = 35-1; %-*-
                Nr_hemi = 35;
            case 'rh'
                %Nr_hemi = 35-1;
                Nr_hemi = 35;
            otherwise
                Nr_hemi = [];
        end
        %lbls = [1:3, 5:35];
        % -*- minus 1 because label 4 (corpuscallosum) not included in
        % cortical surface. Update: actually it is since the surface is
        % closed to form a topology so CC is included. 
        lbls = 1:35;

    case  'Braak5'

        % Braak has six stages.
        %
        % Based on the list of anatomical regions given in Cho2016 for the
        % six stages, and given the labels of the DT atlas (aparc) given by
        % FreeSurfer, I split the DK regions into 6 labels as:
        %
        % [Braak5 label: 1] Braak I-II: not possible to split I and II
        % since they correspond to a single DK region per hemisphere; in
        % fact, Braak I is sub-cortical, I think.
        %
        % [Braak5 label: 2] Braak III
        % [Braak5 label: 3] Braak IV
        % [Braak5 label: 4] Braak V
        % [Braak5 label: 5] Braak VI
        % --[Braak6 label: 6] unassigned regions
        %
        % Note that there are 3 out of 35 DK regions per hemisphere that
        % can not be exclusively assigned to a single Braak stage thus I
        % assign them to a different label (6) as specified above. These
        % regions are:
        % [1] banksts       (DK label 1001/2001) [falls bw Braak IV & V]
        % [2] pericalcarine (DK label 1021/2021) [falls bw Braak V & VI]
        % [3] temporal pole (DK label 1033/2033) [adjacent to Braak IV & V]
        %
        % These regions are thus not assigned to any of the lables in what
        % we denote here as Braak5. 
        % 
        % Note that [1] and [2] are sulci, the banks of which are split
        % between two different stages, whereas for [3] that is part of
        % the temporal lobe it cannot be assigned to one of the temporal
        % cortices (inferior & middle temporal: Braak IV, superior
        % temporal: Braak V).

        ParcName = [];

        switch hemi
            case 'lh'
                Nr_hemi = 5;
            case 'rh'
                Nr_hemi = 5;
            otherwise
                Nr_hemi = [];
        end
        lbls = 1:5;

    case  'Braak3'
        
        % Braak has six stages.
        %
        % We have 3 levels defined by FreeSurfer's aparc+aseg regions as
        % used by BF-folks:
        % Stages I-II   -> my label: 1
        % Stages III-IV -> my label: 2
        % Stages V-VI   -> my label: 3
        % 
        % We cannot seperate e.g. stages I and II since I is sub-cortical.
        % Maybe stages III, IV, V, & VI can be seperated; need to explore. 

        ParcName = [];

        switch hemi
            case 'lh'
                Nr_hemi = 3;
            case 'rh'
                Nr_hemi = 3;
            otherwise
                Nr_hemi = [];
        end
        lbls = 1:3;

        case  'EBM5'
            % Event-Based Modelling (EBM) has five stages.
            %
            % Each stage is specified as a set of FreeSurfer's aparc+aseg
            % regions as used by BF-folks:
            
            ParcName = [];

            switch hemi
                case 'lh'
                    Nr_hemi = 5;
                case 'rh'
                    Nr_hemi = 5;
                otherwise
                    Nr_hemi = [];
            end
            lbls = 1:5;
end

switch atlas
    case 'PALS_B12_Brodmann'
        param.Nr = 141+173;
        
    otherwise
        param.Nr = 2*Nr_hemi;
end
param.Nr_hemi  = Nr_hemi;
param.ParcName = ParcName;
param.lbls     = lbls;
