function [f_surfrois, parcinfo] = hb_aparcfsavg2subj(f_rois,opts,varargin)
% 1. project parcellation atlas* fr fsaverage to subject; .sh script.
% 2. project surface labels to volume.
% 3. assign label to voxels within ribbon not labelled in 2.
% 4. merge lh and rh to get a single labelled map.
% 5. also build 4D version of 4.
%
% * e.g. Glassers, Schaefer's, etc.
%
%
% dependencies:
%   .FreeSurfer
%   .github/aitchbi/matlab-utils
%
% h behjat

p = inputParser;
addParameter(p, 'JustGetSurfaceParcellation', false);
addParameter(p, 'AlsoGetRefRibbon', []); % [*]
addParameter(p, 'JustGetSurfaceParcellationNames', false);
addParameter(p, 'DeleteCopiedFsaverageHB', true);
parse(p,varargin{:});
opts2 = p.Results;
% [*] only applicable if JustGetSurfaceParcellation is true. 

parcinfo = struct;

opts = mergestructs(opts,opts2);

if not(isfield(opts,'OverWriteExistingRois'))
    opts.OverWriteExistingRois = false;
end

if not(isfield(opts,'MinimalOutput'))
    opts.MinimalOutput = true;
end

[RunParcellation, names, Nr_hemi, f_surfrois] = checkstuff(f_rois, opts);

if opts.JustGetSurfaceParcellationNames
    return;
end

if RunParcellation

    if ~isfield(opts,'RunPar')
        opts.RunPar = false;
    end

    if isfield(opts,'DoPruning')
        DoPruning = opts.DoPruning;
    else
        DoPruning = true;
    end

    for hemi = {'lh','rh'}

        switch hemi{:}
            case 'lh'
                n_src = names.n_src_lh;
                f_src = names.f_src_lh;
            case 'rh'
                n_src = names.n_src_rh;
                f_src = names.f_src_rh;
        end

        f_annot = f_surfrois.(hemi{:});

        if opts.JustGetSurfaceParcellation
            if isempty(opts.AlsoGetRefRibbon)
                AGRF = false;
            else
                AGRF = opts.AlsoGetRefRibbon;
            end
            if AGRF
                d1 = exist(f_annot, 'file');
                d2 = exist(names.f_rib_lhrh, 'file');
                if not(d2)
                    if exist([names.f_rib_lhrh, '.gz'], 'file')
                        d2 = true;
                        RefRibNeeded = false;
                    else
                        fprintf('\n..Reference ribbon missing: %s', names.f_rib_lhrh);
                        RefRibNeeded = true;
                    end
                else
                    RefRibNeeded = false;
                end
                chk = d1 && d2;
                t = '(& reference ribbon) exists';
            else
                chk = exist(f_annot, 'file');
                t = 'exists';
                RefRibNeeded = false;
            end

            if chk
                fprintf( ...
                    '\n\n..Surface %s parcellation %s: %s \n', ...
                    hemi{:}, ...
                    t, ...
                    f_annot);
                continue;
            end

        else

            [~,~] = mkdir(fileparts(f_src));

            f_src_proc1 = strrep(f_src,'.nii','_proc1.nii');
            f_src_proc2 = strrep(f_src,'.nii','_proc2.nii');
            f_src_proc3 = strrep(f_src,'.nii','_proc3.nii');
            f_src_proc4 = strrep(f_src,'.nii','_proc4.nii');
        end

        %-Initiate labels--------------------------------------------------

        % copy fsaverage_hb folder to G.f.T1w
        switch hemi{:}
            case 'lh'
                d_fsavg_tmp = fullfile(opts.dirs.subjs,'fsaverage_hb');
                if ~exist(d_fsavg_tmp, 'dir')
                    copyfile(opts.dirs.fsavg_hb, d_fsavg_tmp);
                end
        end

        % project fsaverage atlas to subject
        fprintf(...
            sprintf('\n..Projecting %s to subject space..',...
            opts.WhichParcellation));

         cmd = sprintf(...
            '%s -d %s -s %s -h %s -r %s -p %s -j %s',...
            opts.files.sh_aparcfsavg2subj,...
            opts.dirs.subjs,...
            opts.ID,...
            hemi{:},...
            opts.dirs.freesurfer,...
            names.ParcName, ...
            num2str(opts.JustGetSurfaceParcellation));

        runcmd(cmd,'Error in aparcfsavg2subj.sh.');

        d_niftitmp = fullfile(opts.dirs.subjs,opts.ID,'HB',[n_src,'.nii']);

        d_annottmp = fullfile(opts.dirs.subjs,opts.ID,'HB',[n_src,'.annot']);
       
        parcinfo = verifyannot(d_annottmp, hemi{:}, Nr_hemi, parcinfo);
        
        copyfile(d_annottmp,f_annot);

        delete(d_annottmp);
        
        if opts.DeleteCopiedFsaverageHB
            switch hemi{:}
                case 'rh'
                    rmdir(d_fsavg_tmp,'s');
            end
        end

        if not(opts.JustGetSurfaceParcellation)
            % transfer result
            if exist(d_niftitmp,'file')
                copyfile(d_niftitmp,f_src); % non-ziped
                delete(d_niftitmp);
            else
                dgz = [d_niftitmp,'.gz'];
                assert(logical(exist(dgz,'file')));
                copyfile(dgz,[f_src,'.gz']); % zipped
                delete(dgz);
            end
        end

        %-Prepare ribbon---------------------------------------------------
        f_rib = strrep(f_src,'.nii',...
            '.ref_ribbon_extracted_from_mri_ribbon_mgz.nii');

        fprintf('\n..Preparing %s ribbon..',hemi{:});

        f_ribnii = opts.files.ribbonT1w;

        if not(exist(f_ribnii, 'file'))

            d = strrep(f_ribnii, '.nii', '.mgz');

            f_ribmgz = strrep(d, '.gz', '');

            assert(exist(f_ribmgz, 'file'));

            cmd = sprintf('%s -d %s -r %s -i %s -o %s',...
                fullfile(opts.dir_hbfssh, 'hb_fs_mri_convert.sh'),...
                opts.dirs.subjs,...
                opts.dirs.freesurfer,...
                f_ribmgz,...
                f_ribnii...
                );

            runcmd(cmd,'mgz to nii file conversion failed.');
        end

        [v_ribs,h_ribs] = hb_nii_load(opts.files.ribbonT1w);

        h = struct;
        h.fname = f_rib;
        h.dim = h_ribs.dim;
        h.mat = h_ribs.mat;
        h.dt = [2 0];

        switch hemi{:}
            case 'lh'
                k = 3;
            case 'rh'
                k = 42;
        end

        v = zeros(h.dim);

        v(v_ribs==k)=1;

        spm_write_vol(h,v);

        if opts.JustGetSurfaceParcellation 
            % f_rib & f_src match voxel to voxel [*] so no need for this
            % reslicing, since we also don't have f_src when
            % JustGetSurfaceParcellation=true
            %
            % [*] checked on: 10.05.2024]
        else
            % reslice (downsample) to FreeSurfer's ribbon to match f_src
            hb_nii_reslice(...
                f_rib,...
                f_src,...
                0,...
                f_rib,...
                true);
            
            [v,h] = hb_nii_load(f_rib);
        end

        % make connected
        v = hb_make_connected(v,26);

        spm_write_vol(h,v); % updated f_rib needed for hb_prune_graph.m

        if DoPruning
            % remove voxels that are connected to other voxels only via
            % spurious, anatomically unjustifiable edges
            indV = find(v);

            A = hb_get_adjacency(v,26,'weight','no','verbose',0);

            % prune A
            op = hb_prune_graph(...
                A,...
                f_rib,...
                opts.files.srfs.(hemi{:}),...
                'parallelize',opts.RunPar);

            % set voxels removed during pruning to 0
            d = op.pruning.ind_pre_pruning_A_remained_post_pruning;

            v(indV(not(d)))=0;

            % verify connectedness & write
            v0 = v;

            v = hb_make_connected(v0,26);
            if not(isequal(v,v0))
                warning('As I see, should already be connected, no?');
            end

            spm_write_vol(h,v);
        end

        switch hemi{:}
            case 'lh'
                f_rib_lh = f_rib;
            case 'rh'
                f_rib_rh = f_rib;
        end

        if opts.JustGetSurfaceParcellation
            continue;
        end

        %-Process label----------------------------------------------------
        fprintf('\n..Processing labels..');

        %-Remove labels outside ribbon.
        f_in = f_src;

        f_out = f_src_proc1;

        hb_gunzip(f_src);

        rmnonrib(f_in,f_out,f_rib,Nr_hemi);

        fprintf('\n..Number of labels removed: %d/%d [outside ribbon]',...
            nnz(hb_nii_load(f_in))-nnz(hb_nii_load(f_out)),...
            nnz(hb_nii_load(f_in)));

        %-Fill in holes.
        f_in  = f_src_proc1;

        f_out = f_src_proc2;

        [v,h] = hb_nii_load(f_in);

        v_orig = v;

        v0 = false(size(v));

        c1 = zeros(Nr_hemi,1);

        c2 = zeros(Nr_hemi,1);

        for iR=1:Nr_hemi
            d = v0;
            d(v==iR) = true;
            c1(iR) = nnz(d);
            d = and(imfill(d,'holes'),not(d)); % filled voxels
            d = and(d,not(v_orig)); % only those unlabelled
            c2(iR) = nnz(d);
            v(d) = iR; % assign labels
            showprg(iR,Nr_hemi,0,'..Filling holes.. ');
            %hb_progress(iR,Nr_hemi,'tag','..Filling holes.. ');
        end

        d = nnz(v)-nnz(v_orig);

        fprintf([...
            '\n..Number of voxels filled: %d ',...
            '[%.02f %% increase in lablled voxels]'...
            ],...
            d,d/nnz(v_orig)*100);

        h.fname = f_out;

        spm_write_vol(h,v);

        fprintf(...
            '\n..Number of labels generated: %d [filled holes]',...
            nnz(hb_nii_load(f_out))-nnz(hb_nii_load(f_in)));

        %-Remove labels outside ribbon, again.
        f_in  = f_src_proc2;

        f_out = f_src_proc3;

        rmnonrib(f_in,f_out,f_rib,Nr_hemi);

        fprintf(...
            '\n..Number of labels removed: %d/%d [outside ribbon]',...
            nnz(hb_nii_load(f_in))-nnz(hb_nii_load(f_out)),...
            nnz(hb_nii_load(f_in)));

        %-Label unlabelled voxels within ribbon.
        %
        % 1. For each hemi, build 26-conn adjacency matrix (A) from ribbon.
        % 2. Prune A (uisng pial and white) to just keep anatomically
        % justified edges; i.e., build ribbon voxBG.
        % 2. Store list initially labeled voxels (ILV) and their labels.
        % 3. For each unlabelled voxel (UV) in ribbon:
        %   - extract its neighb voxels (NV) from A.
        %   - extract UL=unique(NV labels), M=length(UL)
        %   - select dominat label based on count.
        % Note that some voxels in ribbon remain unlabled since they
        % are essentially too far away for ILV; these are region in
        % Freesurfer ribbon that are not included in Glasser's atlas,
        % so it makes sense to keep them unlabelled and not use them
        % for defining ROIs when building seed-based FC maps.
        %
        %
        % NOTE: step 2. is optional; only done if DoPruning=true.

        f_in  = f_src_proc3;

        f_out = f_src_proc4;

        rib_lbl = assignlabels(...
            f_in,...
            f_rib,...
            opts.files.srfs.(hemi{:}),...
            opts.RunPar,...
            DoPruning);

        v_rb = hb_nii_load(f_rib);

        v = zeros(size(v_rb));

        v(logical(v_rb)) = rib_lbl;

        [~,h] = hb_nii_load(f_in);

        h.fname = f_out;

        spm_write_vol(h,v);

        fprintf('\n..Number of labels generated: %d %s',...
            nnz(hb_nii_load(f_out))-nnz(hb_nii_load(f_in)),...
            '[based on adjacencies as defined by ribbon voxBG]');

        d = nnz(hb_nii_load(f_rib));

        d1 = d-nnz(hb_nii_load(f_out));

        d2 = d1/d*100;

        fprintf(...
            '\n..Number of voxels in ribbon left unlabelled: %d %s',...
            d1,...
            sprintf('[%.1f%% of total number of voxels in ribbon]',d2));
        if d2>2
            warning(...
                'notable chunk of %s ribbon left unlabelled: %.1f%%', hemi{:}, d2);
        else
            fprintf(...
                '\n..A small portion of %s ribbon left unlabelled: %.1f%% \n', hemi{:}, d2);
        end

        if not(opts.MinimalOutput)
            gzip(f_src);
            gzip(f_src_proc1);
            gzip(f_src_proc2);
            gzip(f_src_proc3);
        end
        delete(f_src);
        delete(f_src_proc1);
        delete(f_src_proc2);
        delete(f_src_proc3);
    end

    if opts.JustGetSurfaceParcellation && not(RefRibNeeded)
        return;
    end

    % Build reference both-hemisphere ribbon ------------------------------
    % f_src/f_rois/opts.files.ribbonT1w have matching dim and mat
    d = spm_vol(opts.files.ribbonT1w);
    info = struct;
    info.dir = names.d_save;
    info.dim = d.dim;
    info.mat = d.mat;
    sts = buildrefrib(info, names.f_rib_lhrh, f_rib_lh, f_rib_rh, [], []);

    delete(f_rib_lh);
    delete(f_rib_rh);

    assert(sts==1);

    if opts.JustGetSurfaceParcellation
        return;
    end

    % Merge lh & rh labels-------------------------------------------------

    fprintf('\n..Merging lh and rh labels..');

    f1 = strrep(names.f_src_lh,'.nii','_proc4.nii');

    f2 = strrep(names.f_src_rh,'.nii','_proc4.nii');

    [v1,h1] = hb_nii_load(f1);

    [v2,h2] = hb_nii_load(f2);

    assert(isequal(h1.dim,h2.dim));
    assert(all(h1.mat(:)-h2.mat(:)< 1e-3));

    I = find(v2);
    v2(I) = v2(I)+ Nr_hemi; %[*] 
    %[*] +Nr_hemi for rh label: 
    % lh label/region R equivalent to rh label/region R+Nr_hemi
    v = v1+v2;
    h_rois = h1;
    h_rois.fname = f_rois;
    spm_write_vol(h_rois,v);

    if not(opts.MinimalOutput)
        gzip(f1);
        gzip(f2);
    end
    delete(f1);
    delete(f2);
end

% Reslice parcellation? ---------------------------------------------------
if isfield(opts.files,'refnii')

    f_refnii = opts.files.refnii;

    if not(isempty(f_refnii))

        if isfield(opts.files,'rois_ref') && not(isempty(opts.files.rois_ref))
            f_rois_ref = opts.files.rois_ref;
        else
            f_rois_ref = strrep(f_rois,'.nii',[opts.resTag,'.nii']);
        end

        [~,n1]  = fileparts(f_rois);
        [p2,n2] = fileparts(f_rois_ref);
        if contains(n2, n1)
            if strfind(n2,n1)==1
                % n2 format: n1 followed by a tag e.g. n1.res1250
                % thus append refnii befor the tag
                n3 = [n1, '_refrib', n2(length(n1)+1:end)];
                f_rib_lhrh_ref = fullfile(p2, [n3, '.nii']);
            else
                error('extend');
            end
        else
            error('extend');
        end

        d1 = not(exist(f_rois_ref,'file'));

        d2 = not(exist([f_rois_ref,'.gz'],'file'));

        d3 = not(exist(f_rib_lhrh_ref,'file'));

        d4 = not(exist([f_rib_lhrh_ref,'.gz'],'file'));

        if and(d1,d2) || and(d3,d4)

            %-Reslice to cerebrum voxBG mask dim (res1250, res2000, etc).
            hb_nii_reslice(...
                f_rois,...
                f_refnii,...
                0,...
                f_rois_ref,...
                true);
            gzip(f_rois_ref);
            delete(f_rois_ref);

            %-Also the asociated refrib.
            hb_nii_reslice(...
                names.f_rib_lhrh,...
                f_refnii,...
                0,...
                f_rib_lhrh_ref,...
                true);
            gzip(f_rib_lhrh_ref);
            delete(f_rib_lhrh_ref);
        end
    end
end

if exist(f_rois, 'file')
    gzip(f_rois);
    delete(f_rois);
end
end

%==========================================================================
function rmnonrib(f_in, f_out, f_rib, Nr_hemi)
[v,h] = hb_nii_load(f_in);
[v_rib,h_rib] = hb_nii_load(f_rib);
verifylabelset(v, Nr_hemi);
assert(isequal(h.dim,h_rib.dim));
assert(all(h.mat(:)-h_rib.mat(:)< 1e-3));
v_bkup = v;
v(not(v_rib)) = 0;
lbl_missing = verifylabelset(v, Nr_hemi);
if not(isempty(lbl_missing))
    %---Generate missing label from neighoring voxels outside ribbon---
    % For each missing label:
    % [1] dilate outside ribbon label
    % [2] check whetehr it bleeds into ribbon
    % [3] repeat [1]-[2] until [2] happens

    se = strel('sphere', 1); % sphere with a radius of one voxel
    for k=1:length(lbl_missing)
        lbl = lbl_missing(k);
        v0 = zeros(size(v));
        v0(v_bkup==lbl) = 1; % single-label map
        DilateMore = true;
        RoundMax = 3;
        Rounds = 1;
        while DilateMore
            if Rounds>RoundMax
                errmsg = ['label number ', ...
                    sprintf('%d (of %d hemi labels)',lbl, Nr_hemi),...
                    'did not bleed into ribbon even after ', ...
                    sprintf('%d rounds of dilation', RoundMax)];
                error('fishy: %s', errmsg);
                % intuitively, 3 dilation rounds should be sufficient,
                % otherwise fishy label position.
            else
                Rounds = Rounds +1;
            end
            v0 = imdilate(v0,se); % dilate
            if nnz(v0(v_rib==1)) % bleeded into ribbon?
                v0(v_rib) = -1;
                v_bkup(v0==-1) = lbl; % NOTE1
                DilateMore = false;
            end
        end
    end
    verifylabelset(v_bkup, Nr_hemi); % NOTE2

    % NOTE1: we might overwrite some other labels, but that should be fine
    % since it's so few voxels that will penerate ribbon with the samll se
    % we use.

    % NOTE2: only way this could potentialy still not be verified is that
    % within the above fro loop some labels of later iteration overwrite
    % those generated in earlier iterations, which would require those
    % missing labels to be exactly adjacent; a quite unlikely scenario.
    v = v_bkup;
    v(not(v_rib)) = 0;
    verifylabelset(v, Nr_hemi);
    disp('');
    disp('all good.');
end
h.fname = f_out;
spm_write_vol(h,v);
end

%==========================================================================
function lbl_missing = verifylabelset(v, Nr_hemi)
if not(isempty(Nr_hemi))
    d = Nr_hemi+1;
    if length(unique(v))~=d
        disp(' ')
        disp('number of unique labels:')
        disp(length(unique(v)))
        disp(' ')
        disp('expected number of unique labels:')
        disp(d)
        disp(' ')
        disp('missing labels:')
        lbl_missing = setdiff(0:Nr_hemi, unique(v))
        Nm = length(lbl_missing);
        if nargout==0
            error(...
                'fishy: %d/%d labels missing (not even outside ribbon)',...
                Nm, Nr_hemi);
        end
    else
        lbl_missing = [];
    end
else
    lbl_missing = [];
end
end

%==========================================================================
function lb_rb = assignlabels(f_lbl,f_rib,f_srfs,RunPar,DoPruning)

[v_rb,h_rb] = hb_nii_load(f_rib);
[v_lb,h_lb] = hb_nii_load(f_lbl);
assert(isequal(h_lb.dim,h_rb.dim));
assert(all(h_lb.mat(:)-h_rb.mat(:)< 1e-3));

I_rb = find(v_rb);
I_lb = find(v_lb);
assert(all(ismember(I_lb,I_rb)),...
    'v_lb should have been masked by v_rb early in the pipeline.');

N_rb = length(I_rb);

I_done = ismember(I_rb,I_lb); % which voxels are already labelled?

lb_rb = zeros(N_rb,1);
lb_rb(I_done) = v_lb(I_lb);

A0 = hb_get_adjacency(v_rb,26,'verbose',0);

assert(isequal(size(A0,1),N_rb),...
    'v_rb should have been made connected earlier in the pipeline.');

if DoPruning

    d = hb_prune_graph(...
        A0,...
        f_rib,...
        f_srfs,...
        'parallelize',RunPar);
    A = d.A;

    if not(isequal(size(A0,1),size(A,1)))
        % set f_rib voxels removed during pruning to 0
        d = d.pruning.ind_pre_pruning_A_remained_post_pruning;
        v_rb(I_rb(not(d))) = 0;
        if 0 % for dubug
            sl = 100;
            hm_plot_adjacency_diff(A0(d,d),A,v_rb,sl);
        end
        spm_write_vol(h_rb,v_rb);

        % update stuff
        I_rb = find(v_rb);
        N_rb = length(I_rb);
        I_done = ismember(I_rb,I_lb); % which voxels are already labelled?
        lb_rb = zeros(N_rb,1);

        whos I_done
        whos lb_lb
        whos v_lb
        whos I_lb

        lb_rb(I_done) = v_lb(I_lb);
    else
        % no update needed on on f_rib
    end
else
    A = A0;
end

assert(isequal(size(A,1),length(I_done)),'fishy: size mismatch');

E = A(not(I_done),:); % Neighboring Voxels of Unlabelled Voxels (NVoUV)
L = E.*(lb_rb');      % labels of NVoUV

clear A0 A E

N_tlb = nnz(not(I_done));

if RunPar
    [y,p] = hb_spmd_prepare({1:N_tlb},1);
    spmd
        d = selectlabels(L,y{1}{labindex},0);
    end
    lbls = zeros(N_tlb,1);
    for k=1:p.NumWorkers
        lbls(y{1}{k}) = d{k};
    end
    delete(p);
else
    lbls = selectlabels(L,1:N_tlb,1);
end
lb_rb(not(I_done)) = lbls;

end

%==========================================================================
function lbls = selectlabels(L,rows,ShowPrgs)
N = length(rows);
lbls = zeros(N,1);
for iR=1:N
    r = rows(iR);
    l = L(r,logical(L(r,:)));
    u = unique(full(l));
    u(u==0) = []; % candidate labels
    if isempty(u)
        continue;
    end
    c = zeros(1,length(u));
    for iU=1:length(u)
        c(iU) = nnz(l==u(iU));
    end
    [c,d] = sort(c,'descend'); %#ok<ASGLU>
    lbls(iR) = u(d(1)); % the dominant candidate label
    if ShowPrgs
        if or(or(iR==1,iR==N), rem(iR,100)==0)
            showprg(iR,N,0,'..Assigning label to unlabelled voxels.. ');
        end
    end
end
end

%==========================================================================
function runcmd(cmd,errmsg)
[sts,log] = system(cmd);
if sts==0
    return;
end
sprintf('*** system run command log: \n\n');
log %#ok<NOPRT>
error(errmsg)
end

%==========================================================================
function showprg(n,N,init,tag)
l = numel(num2str(N));
if n==1 || init
    fprintf(['\n',tag]);
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end

%==========================================================================
function sts = buildrefrib(info, f_rib_lhrh, f_lh, f_rh, WhichAtlas, atlases)
% a temp fix, to maybe remove later; this is to prevent rebuilding
% parcellations in case they exist but only the ref ribbon is missing since
% it has not been generated from existing f_lh and f_rh.
if isempty(f_lh) || isempty(f_rh)

    t2 = 'ref_ribbon_extracted_from_mri_ribbon_mgz';

    switch WhichAtlas
        case atlases
            d = WhichAtlas(9:11);
            t1 = sprintf('Schaefer2018_%sParcels_7Networks_order', d);
            f_lh = fullfile(info.dir, sprintf('lh.%s.%s.nii.gz',t1,t2));
            f_rh = fullfile(info.dir, sprintf('rh.%s.%s.nii.gz',t1,t2));
        otherwise
            error('extend');
    end
end

d1 = not(exist(f_lh,'file'));
d2 = not(exist(f_rh,'file'));
if d1 || d2
    d1 = exist([f_lh, '.gz'],'file');
    d2 = exist([f_rh, '.gz'],'file');
    if d1 && d2
        gunzip([f_lh, '.gz']);
        gunzip([f_rh, '.gz']);
    else
        sts = 0;
        return;
    end
end

h = struct;
h.fname = f_rib_lhrh;
h.dim = info.dim;
h.mat = info.mat;
h.dt = [2 0];
h = spm_create_vol(h);

v = zeros(h.dim);

h_lh = spm_vol(f_lh);

h_rh = spm_vol(f_rh);

assert(isequal(h_lh.dim,h_rh.dim));

assert(all((h_lh.mat(:)-h_rh.mat(:))<1e-3));

v_lh = logical(spm_read_vols(h_lh));

v_rh = logical(spm_read_vols(h_rh));

assert(isempty(intersect(find(v_lh),find(v_rh))));

v(v_lh) = 3;

v(v_rh) = 42;

spm_write_vol(h,v);

gzip(f_rib_lhrh);

delete(f_rib_lhrh);

sts = 1;

end

%==========================================================================
function opts = mergestructs(opts1, opts2)
F1 = fieldnames(opts1);
F2 = fieldnames(opts2);
assert(~any(ismember(F1, F2)), 'Fishy function call.');
opts = struct;
for k=1:length(F1)
    f = F1{k};
    opts.(f) = opts1.(f);
end
for k=1:length(F2)
    f = F2{k};
    opts.(f) = opts2.(f);
end
end

%==========================================================================
function [RunParc, names, Nr_hemi, f_surfrois] = checkstuff(f_rois, opts)

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

if contains(f_rois, '.gz')
    f_roisgz = f_rois;
    f_rois = strrep(f_rois, '.gz', '');
elseif endsWith(f_rois, '.nii')
    f_roisgz = strrep(f_rois, '.nii', '.nii.gz');
end

names = struct;

[p,n] = fileparts(f_rois);

names.f_rib_lhrh = fullfile(p,[n,'_refrib.nii']);
f_rib_lhrh_gz    = strcat(names.f_rib_lhrh,'.gz');

if opts.JustGetSurfaceParcellation

    RunParc = true;

else

    if exist(f_roisgz,'file')
        if opts.OverWriteExistingRois
            delete(f_roisgz);
            FileExists = false;
        else
            FileExists = true;
        end
    elseif exist(f_rois,'file')
        if opts.OverWriteExistingRois
            delete(f_rois);
            FileExists = false;
        else
            FileExists = true;
            gzip(f_rois);
            delete(f_rois);
        end
    else
        FileExists = false;
    end

    if FileExists
        if not(exist(f_rib_lhrh_gz, 'file'))
            d = spm_vol(f_rois);
            infofrois = struct;
            infofrois.dir = fileparts(f_rois);
            infofrois.dim = d.dim;
            infofrois.mat = d.mat;
            sts = buildrefrib(infofrois, names.f_rib_lhrh, [], [], opts.WhichParcellation, SchaeferYeo7);
            if sts==0
                fprintf('\n\n..Parcellation exists but ref ribbon mising: %s \n', names.f_rib_lhrh);
                fprintf('\n..Recomuting parcellation.. ');
            end
            RunParc = true;
        else
            fprintf('\n\n..Parcellation exists: %s', f_roisgz);
            RunParc = false;
        end
    else
        RunParc = true;
    end
end

SchaeferYeo = [SchaeferYeo7; SchaeferYeo17];

getSchYeo7 = @(x) sprintf('Schaefer2018_%dParcels_7Networks_order',x);

getSchYeo17 = @(x) sprintf('Schaefer2018_%dParcels_17Networks_order',x);

switch opts.WhichParcellation

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
        d_end = strfind(opts.WhichParcellation,'Yeo')-1;
        Nr_hemi = str2double(opts.WhichParcellation(d_str:d_end))/2;

        d = Nr_hemi*2;

        switch opts.WhichParcellation
            case SchaeferYeo7
                ParcName = getSchYeo7(d);

            case SchaeferYeo17
                ParcName = getSchYeo17(d);
        end

    otherwise
        error('extend code.');
end

d_save = fileparts(f_rois);

[~,~] = mkdir(d_save);

names.ParcName = ParcName;
names.n_src_lh = sprintf('lh.%s',ParcName);
names.n_src_rh = sprintf('rh.%s',ParcName);
names.f_src_lh = fullfile(d_save,[names.n_src_lh,'.nii']);
names.f_src_rh = fullfile(d_save,[names.n_src_rh,'.nii']);
names.d_save = d_save;

f_surfrois = struct; % with fields lh & rh
for hemi = {'lh','rh'}
    [p,n,e] = fileparts(f_rois);
    assert(isequal(e,'.nii'));
    f_surfrois.(hemi{:}) = fullfile(p,[hemi{:}, '.', n, '.annot']);
end
end

%==========================================================================
function parcinfo = verifyannot(f, hemi, N, parcinfo)
assert(exist(f, 'file'), 'annotation file missing');
[~, s, tbl] = read_annotation(f);
L = unique(s);
M = tbl.numEntries;
U = length(L);
assert(isequal(U, M), 'fishy');
assert(isequal(N, M-1), 'fishy'); % one label is background 
I_bg = 1; % first label in table is background
L_bg = 65793; % background label; RGB [1 1 1]
assert(tbl.table(I_bg,5)==L_bg, 'Background label not found');
assert(contains(tbl.struct_names{I_bg}, 'Background'));
assert(contains(tbl.struct_names{I_bg}, 'Medial_Wall'));
C_bg  = nnz(s==L_bg);
n_bg  = tbl.struct_names{I_bg};
C_lbl = zeros(N,1);
n_lbl = cell(N,1);
for k=2:M
    j           = k-1;
    l           = tbl.table(k,5);
    C_lbl(j)    = nnz(s==l);
    n_lbl{j} = tbl.struct_names{j};
end
n = nnz(C_lbl==0);
if n~=0
    error('Subject-specific parcels missing: %d/%d', n, N);
end
[Nmin, Imin] = min(C_lbl);
[Nmax, Imax] = max(C_lbl);
lblmin = n_lbl{Imin};
lblmax = n_lbl{Imax};
t = 'Number of surface vertices';

disp(' ');
fprintf('\n.Subject-specific %s surface parcellation done.', hemi);
fprintf('\n..%s in smallest parcel [%s]: %d (parcel #%d: %s)', t, hemi, Nmin, Imin, lblmin);
fprintf('\n..%s in largest parcel  [%s]: %d (parcel #%d: %s)', t, hemi, Nmax, Imax, lblmax);
fprintf('\n..%s in parcels [%s, median]: %d', t, hemi, round(median(C_lbl)));
disp(' ');

parcinfo.(hemi).labels_num  = 1:N;
parcinfo.(hemi).labels_name = n_lbl;
parcinfo.(hemi).labels_size = C_lbl;
parcinfo.(hemi).labels_cmap = tbl.table(2:end,5);
parcinfo.(hemi).bg_name     = n_bg;
parcinfo.(hemi).bg_size     = C_bg;
parcinfo.(hemi).bg_cmap     = tbl.table(1,5);
end