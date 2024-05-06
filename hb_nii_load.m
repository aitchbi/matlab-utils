function [v,h] = hb_nii_load(f,varargin)
%HB_NII_LOAD load nifti volume and header using external software. 

d = inputParser;
addParameter(d,'JustGetHeader', false);
addParameter(d,'IndicesToLoad', []);
addParameter(d,'FramesToLoad', []);
addParameter(d,'HeaderType', 'spm');
addParameter(d,'LoadAsVectors', false);
addParameter(d,'DuplicateThenUnzip', false);
parse(d,varargin{:});
opts = d.Results;

assert(ischar(f),'Enter absolute address of nifti file.');

if contains(f,'.gz')
    if exist(f,'file')
        fgz = f;
        f = dounzip(fgz, opts.DuplicateThenUnzip);
        DELNONZIP = true;
    elseif exist(strrep(f,'.gz',''),'file')
        DELNONZIP = false;
    else
        error('File missing: %s',f);
    end
   
else
    if exist(f,'file')
        DELNONZIP = false;
    else
        fgz = [f,'.gz'];
        if exist(fgz,'file')
            f = dounzip(fgz, opts.DuplicateThenUnzip);
            DELNONZIP = true;
        else
            error('File missing: %s',f);
        end
    end
end

if isempty(opts.IndicesToLoad)
    FASTLOAD = false;
else
    FASTLOAD = true;
end
    
switch opts.HeaderType
    case 'spm'
        %assert(exist('spm_vol.m','file'),'"spm" package not in path.')
        h = spm_vol(f);
        if opts.JustGetHeader
            v = [];
            return;
        end
        Nv = length(h);
        if isempty(opts.FramesToLoad)
            frames = 1:Nv;
        else
            frames = opts.FramesToLoad;
        end
        Nf = length(frames);
        assert(Nf<=Nv);
        assert(all(ismember(frames,1:Nv)));
        if Nv==1
            v = spm_read_vols(h);
        else
            if FASTLOAD
                h1 = h(1);
                I = opts.IndicesToLoad;
                [x,y,z] = ind2sub(h1.dim,I);
                v1d0 = zeros(prod(h1.dim),1);
                if opts.LoadAsVectors
                    v = zeros(length(I),Nf);
                else
                    v = zeros([h1.dim,Nf]);
                end
                for iF=1:Nf
                    iV = frames(iF);
                    v1d = v1d0;
                    v1d(I) = spm_sample_vol(h(iV),x,y,z,0);
                    if opts.LoadAsVectors
                        v(:,iF) = v1d(I);
                    else
                        v(:,:,:,iF) = reshape(v1d,h1.dim);
                    end
                    if Nf>1
                        prgs(iF,Nf);
                    end
                end
            else
                fprintf('\n..Loading 4D nifti.. \n');
                v = spm_read_vols(h);
                if not(isequal(Nf,Nv))
                    v = v(:,:,:,frames);
                end
            end
        end
        if DELNONZIP
            delete(f);
        end
    otherwise
        error('extend.')
end
end

%==========================================================================
function prgs(n,N)
l = numel(num2str(N));
if n==1
    fprintf('\n..Loading 4D nifti.. ');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end

%==========================================================================
function t = getrandtag
t = sprintf('___tmp___%d',round(rand*1e12));
end

%==========================================================================
function f = dounzip(fgz, DuplicateThenUnzip)
if DuplicateThenUnzip
    fgztmp = strrep(fgz, '.nii', sprintf('%s.nii', getrandtag));
    copyfile(fgz,fgztmp);
    fgz = fgztmp;
end
gunzip(fgz);
if DuplicateThenUnzip
    delete(fgztmp);
end
f = strrep(fgz,'.gz','');
end
