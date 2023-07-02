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
        
    case 'Yeo17'
        Nr_hemi = 17;
        ParcName = 'Yeo2011_17Networks_N1000';
        
    case 'Yeo7'
        Nr_hemi = 7;
        ParcName = 'Yeo2011_7Networks_N1000';
        
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
end

switch atlas
    case 'PALS_B12_Brodmann'
        param.Nr = 141+173;
        
    otherwise
        param.Nr = 2*Nr_hemi;
end
param.Nr_hemi = Nr_hemi;
param.ParcName = ParcName;
