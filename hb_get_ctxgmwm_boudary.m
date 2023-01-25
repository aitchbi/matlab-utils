function f = hb_get_ctxgmwm_boudary(d_fsmri,SaveAuxFiles)
% HB_GET_CTXGMWM_BOUNDARY extracts touching gray matter (GM) and white
% matter (WM) voxels at the boundary of cerebral cortex, that is, the most
% superficial WM voxels, and the inner layer voxels of cerebral cortex. The
% function works on FreeSurfer extracted files.
%
% Input:
%
%   d_fsmri: absolute address of FreeSurfer 'mri' folder for a subject.
%
%   SaveAuxFiles: (optional) 3 auxiliary files are save:
%   1. cerebral cortex GM   
%   2. sub-cortical cerebral GM
%   3. cerebral WM
%
% Output:
%
%   f: absolute address of saved nifti volumes,
%   f.gm_wmborder: cerebral cortex GM border with WM
%   f.wm_gmborder: WM border with cerebral cortex
%   f.gm_ctx     : cerebral cortex GM (optionally saved)
%   f.gm_subctx  : sub-cortical cerebral GM (optionally saved)
%   f.wm_cerebrum: cerebral WM (optionally saved)
%
% Requirements:
%
%   SPM12 software, to read/write nifti files.
%   https://www.fil.ion.ucl.ac.uk/spm/software/spm12
% 
% Hamid Behjat

if ~exist('SaveAuxFiles','var')
    SaveAuxFiles = false;
end

f_rb = fullfile(d_fsmri,'ribbon.nii');     

f_aa = fullfile(d_fsmri,'aparc+aseg.nii'); 

%-Check exist.
chk_exist(f_rb,f_aa);

%-ribbon.
h_rb = spm_vol(f_rb);

v_rb = spm_read_vols(h_rb);

%-aparc+aseg.
v_aa = spm_read_vols(spm_vol(f_aa));

%-REF header.
hh = struct();

hh.dim = h_rb.dim;

hh.mat = h_rb.mat;

hh.dt = [spm_type('uint8'),0];

%-Extract cortex GM.
v1 = ismember(v_rb,[3,42]);

if SaveAuxFiles
    h1 = hh;
    
    h1.fname = strrep(f_rb,'ribbon','cerebrum_gm_ctx');
    
    spm_write_vol(h1,v1);
end

%-Extract cerebral WM.
[v3,h3,h2] = get_wm(v_aa,v1,f_rb,hh,SaveAuxFiles);

%-Strutural element for defining border voxels.
d = ones(3,3,3);

%-GM border with WM.
v4 = and(v1,imdilate(v3,d));

h4 = hh;

h4.fname = strrep(f_rb,'ribbon','cerebrum_gm_wmborder');

spm_write_vol(h4,v4);

%-WM border with GM.
v5 = and(v3,imdilate(v1,d));

h5 = hh;

h5.fname = strrep(f_rb,'ribbon','cerebrum_wm_gmborder');

spm_write_vol(h5,v5);

%-Files names.
f = struct;

f.gm_wmborder = h4.fname; % cerebral cortex GM border with WM

f.wm_gmborder = h5.fname; % WM border with cerebral cortex

if SaveAuxFiles
    f.gm_ctx      = h1.fname; % cerebral cortex GM

    f.gm_subctx   = h2.fname; % sub-cortical cerebral GM
    
    f.wm_cerebrum = h3.fname; % cerebral WM
end
end

%==========================================================================
function chk_exist(f_rb,f_aa)

t = sprintf([
    'is missing. \n. To get it, in Terminal, ',...
    'cd to subject''s FreeSurfer''s mri folder. \n',...
    '. Then run: mri_convert'
    ]);

if ~exist(f_rb,'file')
    
    assert(exist(strrep(f_rb,'.nii','.mgz'),'file'), 'ribbon.mgz missing');
    
    d = 'ribbon';
    
    fprintf('\n. File: %s.nii %s %s.mgz %s.nii \n--------------\n',d,t,d,d);
    
    error('ribbon.nii missing');
end

if ~exist(f_aa,'file')
    
    assert(exist(strrep(f_aa,'.nii','.mgz'),'file'), 'ribbon.mgz missing');
    
    d = 'aparc+aseg';
    
    fprintf('\n. File: %s.nii %s %s.mgz %s.nii \n--------------\n',d,t,d,d);
    
    error('aparc+aseg.nii');
end

if ~exist('spm_vol.m','file')
    
    fprintf('\n. SPM12 software not in path.\n');
    
    fprintf([
        '. Addpath if available, otherwise, first download from here: ',...
        'https://www.fil.ion.ucl.ac.uk/spm/software/spm12\n-------------\n'
        ]);
    
     error('SPM12 software not in path.');
end
end

%==========================================================================
function [v3,h3,h2] = get_wm(v_aa,v1,f_rb,hh,SaveAuxFiles)

%-GM sub cortex.
d1 = [9:13,17:20,26:28, 48:56, 58:60];       % sub-cortical cerebral GM

d2 = [80,81,82];                             % non-WM hypointensities

v2 = ismember(v_aa,[d1,d2]);

d1 = [2,41, 85, 250, 251:255];               % Cerebrum WM; see NOTE3.

d2 = [77,78,79];                             % WM hypointensities

v = ismember(v_aa,[d1,d2]);

d1 = [4,5,43,44,14,15,72, 24, 31,63, 30,62]; % ventricles, CSF, choroid-plexus & vessels

d2 = [1,6,40,45,21:23,25,57,75,76];          % irrelevant

d3 = [7,8,46,47, 16];                        % cerebellum & brain stem; see NOTE2.

v_omit_p1 = ismember(v_aa,[d1,d2,d3]);

v = and(v,not(v_omit_p1));

v = and(v,not(v1));

v = and(v,not(v2));

% Handle unassigned voxels; see NOTE1.
d = [1000:1035, 2000:2035];

v_ua = and(ismember(v_aa,d),not(v1));

v = or(v,v_ua);

v3 = make_connected(v);

if SaveAuxFiles
    h3 = hh;
    
    h3.fname = strrep(f_rb,'ribbon','cerebrum_wm');
    
    spm_write_vol(h3,v3);
    
    h2 = hh;
    
    h2.fname = strrep(f_rb,'ribbon','cerebrum_gm_subctx');
    
    spm_write_vol(h2,v2);
else
    h2 = [];
    h3 = [];
end
end

%==========================================================================
function v = make_connected(v)
% Make array connected. Largest component kept, others deleted.

CC = bwconncomp(v,26);

[~,imax] = max(cellfun(@length,CC.PixelIdxList));

v = zeros(size(v));

v(CC.PixelIdxList{imax}) = 1;
end

%==========================================================================
%
% NOTE1--------------------------------------------------------------------
% There is notable mismatch between the definition of cerebral cortex given
% by FreeSurfer's ribbon and aparc+aseg files.
%
% These are UnAssigned (ua) voxels (quite many) at the gray/white cerebral
% cortex boundary (CLASS1), also, at the pial/CSF boundary (CLASS2), that
% are assigned as being cortex in aparc+aseg.nii, but are not defined as
% cortex in ribbon.nii.
%
% One should note that the cerebral cortex as defined by ribbon is
% more accurate/smooth/neat than that defined by aparc+aseg. So nonw of
% these unassigned voxels will be labelled as cerebral cortex. They will
% either be considered as white matter (CLASS1) or will be dicarded
% (CLASS2).
%
% To handle these voxels, first, all CLASS1 & CLASS2 voxels are assigned as
% being WM. Then, CLASS2 voxels are cleaned out form the WM mask via making
% it connected.
%
% NOTE2--------------------------------------------------------------------
% These should not be in subcortical mask extracted from ribbon, but we
% omit them just in case they are partially there
%
% NOTE3--------------------------------------------------------------------
% Label descriptions:
%
%-Cortical GM: 
%  3  Left-Cerebral-Cortex
% 42  Right-Cerebral-Cortex 
% Note that in aparc+aseg, 3 & 42 are replcaed with values 1xxx & 2xxx.
%
%-Cerebrum WM: 
%   2  Left-Cerebral-White-Matter
%  41  Right-Cerebral-White-Matter
%  85  Optic-Chiasm               [where optic nerves cross, in WM]
% 250  Fornix                     [missing in 100307]
% 251  CC_Posterior
% 252  CC_Mid_Posterior
% 253  CC_Central
% 254  CC_Mid_Anterior
% 255  CC_Anterior 
%
%-Sub-cortical GM:
%  9  Left-Thalamus [missing? at least in 100307]                           
% 10  Left-Thalamus-Proper*                   
% 11  Left-Caudate                            
% 12  Left-Putamen                            
% 13  Left-Pallidum                           
% 17  Left-Hippocampus                        
% 18  Left-Amygdala                           
% 19  Left-Insula           [missing in 100307]                            
% 20  Left-Operculum        [missing in 100307]        
% 26  Left-Accumbens-area                     
% 27  Left-Substancia-Nigra [missing in 100307]                   
% 28  Left-VentralDC 
%
% 48  Right-Thalamus                          
% 49  Right-Thalamus-Proper*                  
% 50  Right-Caudate                           
% 51  Right-Putamen                           
% 52  Right-Pallidum                          
% 53  Right-Hippocampus                       
% 54  Right-Amygdala                          
% 55  Right-Insula           [missing in 100307]                           
% 56  Right-Operculum        [missing in 100307]                         
% 58  Right-Accumbens-area                    
% 59  Right-Substancia-Nigra [missing in 100307]                  
% 60  Right-VentralDC
%
%-Regions to exclude [non-tissue]:
%  4  Left-Lateral-Ventricle
%  5  Left-Inf-Lat-Vent
% 43  Right-Lateral-Ventricle
% 44  Right-Inf-Lat-Vent
% 14  3rd-Ventricle                           
% 15  4th-Ventricle    
% 72  5th-Ventricle 
%
% 24  CSF
%
% 31  Left-choroid-plexus
% 63  Right-choroid-plexus
%
% 30  Left-vessel
% 62  Right-vessel 
%
%-Regions to exclude [tissue, but not included in 'cerebrum' voxBG]:
%  7  Left-Cerebellum-White-Matter
%  8  Left-Cerebellum-Cortex
% 46  Right-Cerebellum-White-Matter
% 47  Right-Cerebellum-Cortex
% 16  Brain-Stem     
%
%-Regions to exclude: 
% [These labels are irrelevant/crap, & typically missing on aseg]
%  1  Left-Cerebral-Exterior    [?; missing in 100307]
%  6  Left-Cerebellum-Exterior  [?; missing in 100307]
% 40  Right-Cerebral-Exterior   [?; missing in 100307]
% 45  Right-Cerebellum-Exterior [?; missing in 100307]
% 21  Line-1                    [?; missing in 100307]
% 22  Line-2                    [?; missing in 100307]
% 23  Line-3                    [?; missing in 100307]
% 25  Left-Lesion               [abnormality; missing in 100307]     
% 57  Right-Lesion              [abnormality; missing in 100307]  
% 75/76                         [FS description: removed. duplicates of 4/43; missing in 100307]
%
%-Regions of unknown nature:
% [they are still tissue, thus, included in WM & subx-ctx GM]
% 77  WM-hypointensities           [in 100307]
% 78  Left-WM-hypointensities      [missing in 100307]
% 79  Right-WM-hypointensities     [missing in 100307]
% 80  non-WM-hypointensities       [in 100307]
% 81  Left-non-WM-hypointensities  [missing in 100307]
% 82  Right-non-WM-hypointensities [missing in 100307]
%
%-Regions that are undetermined: [WHAT TO DO with ???]
% 29  Left-undetermined  [missing in 100307]
% 61  Right-undetermined [missing in 100307]
%
%-Regions that are apparently cerebral cortex, which we will keep: 
% 32  Left-F3orb  [this is seemingly called the pars orbitalis]
% 64  Right-F3orb [this is seemingly called the pars orbitalis]
%
%-Regions, the nature of which is unkown to me:
% [they are typically missing]
% 33  Left-lOg
% 34  Left-aOg
% 35  Left-mOg
% 36  Left-pOg
% 37  Left-Stellate
% 38  Left-Porg
% 39  Left-Aorg
%
% 65  Right-lOg
% 66  Right-aOg
% 67  Right-mOg
% 68  Right-pOg
% 69  Right-Stellate
% 70  Right-Porg
% 71  Right-Aorg
%
% 73  Left-Interior
% 74  Right-Interior
%
% 83  Left-F1
% 84  Right-F1
%
%--------------------------------------------------------------------------
