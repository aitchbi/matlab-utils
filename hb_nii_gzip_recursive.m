function files = hb_nii_gzip_recursive(p,varargin)

d = inputParser;
addParameter(d,'OverwriteExistingGzip',false);
addParameter(d,'JustGetFileList',false); % no gzip done
addParameter(d,'ExcludeEndsWith',[]); % [*]
addParameter(d,'opts',[]); % [**]
addParameter(d,'UnzipGzipniiZip', false); % !!caution!! [***]
parse(d,varargin{:});
opts = d.Results;

opts = fixopts(opts);

% [*] cell array of strings e.g. {'.dscalar.nii', '.dlabel.nii'}
%
% [**] for in-function call see below
%
% [***] if true, recursively unzips existing .zip files, gzips any .nii
% from extracted files, leaves remaining files as they are, zips again

DryRun = opts.JustGetFileList;
OWRITE = opts.OverwriteExistingGzip;

F = dir(p);

for k=1:length(F)
    
    fk = F(k);
    
    f = fullfile(fk.folder, fk.name);

    if endsWith(fk.name, '.nii')
        if exist(f,'file')
            if ~isempty(opts.ExcludeEndsWith)
                SkipFile = 0;
                X = opts.ExcludeEndsWith;
                for i=1:length(X)
                    if endsWith(f, X{i})
                        SkipFile = 1;
                        break;
                    end
                end
                if SkipFile
                    continue;
                end
            end
            if not(DryRun)
                fgz = [f, '.gz'];
                if exist(fgz, 'file')
                    if OWRITE
                        delete(fgz);
                        gzip(f);
                    end
                else
                    gzip(f);
                end
                assert(exist(fgz,'file'));
                delete(f);
            end
            if exist('files', 'var')
                files{length(files)+1} = f; %#ok<AGROW> 
            else
                files = {f};
            end
        end

    elseif fk.isdir
        ok = chkdir(fk.name);
        if ok
            d = hb_nii_gzip_recursive(f, 'opts', opts);
            if not(isempty(d))
                d = d';
                if exist('files', 'var')
                    files = [files, d]; %#ok<AGROW>
                else
                    files = d;
                end
            end
        end

    elseif endsWith(fk.name, '.zip')

        if opts.UnzipGzipniiZip
        
            if exist(f,'file')
                
                % extract zip file into a tmp folder
                d_tmp = fullfile(fk.folder, 'tmpdir_unzip');
                unzip(f, d_tmp);

                % gzip any nifti files
                % note that this is done recursively, so if there is
                % another zip file in the list of extracted files from the
                % above zip file, that will unziped and processed as part
                % of hb_nii_gzip_recursive call
                hb_nii_gzip_recursive(d_tmp, 'opts', opts);

                % delete original zip file 
                delete(f);
                
                % zip the updated contents of the tmp folder
                F_tmp = dir(d_tmp);
                files_tozip = [];
                for k_tmp=1:length(F_tmp)
                    fk_tmp = F(k_tmp);
                    if ismember(fk_tmp.name, {'.', '..', '.DS_Store'})
                        continue;
                    end
                    f_tozip = fullfile(fk_tmp.folder, fk_tmp.name);
                    if exist(f_tozip, 'file')
                        if isempty(files_tozip)
                            files_tozip{1} = f_tozip;
                        else
                            files_tozip{length(files_tozip)+1} = f_tozip; %#ok<AGROW> 
                        end
                    end
                end
                zip(f, files_tozip); % zip files to same name as before

                % delete tmp folder and contents
                for k_dlt=1:length(files_tozip)
                    f_dlt = files_tozip{k_dlt};
                    if exist(f_dlt,'file')
                        delete(f_dlt);
                    end
                end

                % delete tmp folder
                rmdir(d_tmp, 's');
            end
        end
    end
end
if exist('files', 'var')
    files = files';
else
    files = [];
end
end

%==========================================================================
function opts = fixopts(opts)
if isempty(opts.opts)
    return;
else
    F = fieldnames(opts.opts);
    for k=1:length(F)
        f = F{k};
        if not(strcmp(f,'opts'))
            opts.(f) = opts.opts.(f);
        end
    end
    opts.opts = [];
end
end

%==========================================================================
function ok = chkdir(n)
if startsWith(n ,'.')
    ok = 0;
else
    ok = 1;
end
% not going in hidden/system files like:
% '.'
% '..'
% '.DS_Store'
% etc
end