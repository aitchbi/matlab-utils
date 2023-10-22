function hb_nii_write(h,v,varargin)
% HB_NII_WRITE writes nifti file and handles gzip.
%
% Hamid Behjat

d = inputParser;
addParameter(d,'HeaderType', 'spm');
addParameter(d,'DoGzip', []);
addParameter(d,'DeleteAfterGzip', true);
parse(d,varargin{:});
opts = d.Results;

f = h.fname;

if endsWith(f, '.nii.gz')
    DoGzip = true;
    f = strrep(f, '.gz', '');
else
    assert(endsWith(f, '.nii'));
end
h.fname = f;

switch opts.HeaderType
    case 'spm'
        spm_write_vol(h, v);
    otherwise
        error('extend');
end

if isempty(opts.DoGzip)
    %-not specified as optional input.
    %-thus, gzip if filename had .gz. 
    if DoGzip
        dogzip(h.fname, opts.DeleteAfterGzip);
    end
else
    %-gzip status specified as optional input.
    %-disregard filename in header.
    %-i.e, no gzip even if DoGzip=true but opts.DoGzip = false. 
    if opts.DoGzip
        dogzip(h.fname, opts.DeleteAfterGzip);
    end
end

end

%==========================================================================
function dogzip(f, dlt)
gzip(f);
if dlt
    delete(f);
end
end
