function hb_plotoverlay(imgs,plane,slices,ulw,hF,clmap,refimg,putlabels)
% HB_PLOTOVERLAY wrapper function using SPM's slover object, with
% suitable settings, in terms of the colormaps, for overlaying a volume
% with potentially +/- values on a grayscale image.
% 
% Inputs:
%   imgs: cell array of nifti file names. imgs{1} is the overlay image and 
%   is assumed to potentially have both +/- values, like an eigenmode, but 
%   also works very nicely if the image only has + or - values, like the 
%   first Laplacian eigenmode. imgs{2} is the underlay image and it is 
%   assumed to be a grayscale image, i.e., only with + values; for example,
%   a t1w image. 
%   plane: 'axial','sagittal','coronal'
%   slices: vector of slice numbers.
%   ulw: (optional) underlay image weight, scalar in [0,1].
%   hf: (optional) figure handle to plot in. 
%   clmap: (optional) colormap for overlay. 
%
% Dependencies: 
%   SPM12
%
% Hamid Behjat

if ~exist('ulw','var')||isempty('ulw')
    ulw = 0.5;
end
if ~exist('hF','var')||isempty('hF')
    hF = [];
end
if ~exist('clmap','var')||isempty('clmap')
    clmap = [];
end
if ~exist('refimg','var')||isempty(refimg)
    refimg = 1;
end
if ~exist('putlabels','var')||isempty(putlabels)
    putlabels = 1;
end
if ischar(imgs)
    imgs = cellstr(imgs);
end

obj = slover;
obj.cbar = [];

% overlay [eigenvector]
obj.img(1).vol = spm_vol(imgs{1});
obj.img(1).type = 'truecolour';
[mx,mn]  = slover('volmaxmin',obj.img(1).vol);
drange   = [mn,mx];
if isempty(clmap)
    clmap = 'flow.lut';
end
dcmap = clmap;
obj.cbar = [obj.cbar 1];
obj.img(1).cmap  = slover('getcmap',dcmap);
switch clmap
    case 'flow.lut'
        obj.img(1).range = max(abs(drange))*[-1;1];
    otherwise
        obj.img(1).range = drange;
end
obj.img(1).prop  = 1;

% underlay [structural]
if length(imgs)==2
    obj.img(2).vol = spm_vol(imgs{2});
    obj.img(2).type = 'truecolour';
    obj.img(2).cmap = gray;
    [mx,mn] = slover('volmaxmin', obj.img(2).vol);
    obj.img(2).range = [mn mx];
    obj.img(2).prop = ulw;
end
if ~isempty(plane)
    obj.transform = plane;
end
if isempty(hF)
    obj.figure = spm_figure('GetWin', 'Graphics');
else
    obj.figure = hF;
end
obj = fill_defaults(obj);
if refimg~=1
    V = obj.img(refimg).vol;
    D = V.dim(1:3);
    T = obj.transform * V.mat;
    vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
        1 1 D(3); D(1) 1 D(3); 1 D(2:3) ; D(1:3)]';
    corners = T * [vcorners; ones(1,8)];
    SC = sort(corners, 2);
    vxsz = sqrt(sum(T(1:3,1:3).^2));
    obj.slicedef = [SC(1,1) vxsz(1) SC(1,8);SC(2,1) vxsz(2) SC(2,8)];
end
if ~isempty(slices)
    obj.slices = obj.slices(slices);
end
if ~putlabels
    obj.labels = 'none';
end
paint(obj)

