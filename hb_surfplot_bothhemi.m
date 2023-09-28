function opts = hb_surfplot_bothhemi(x,atlas,varargin)
% HB_SURFPLOT_BOTHHEMI same as HB_SURFPLOT.m but showing projection of data
% on both left and right hemipheres in one figure.
%
% *******************************************************
% ***** "preliminary" help; not accurate & thorough *****
% *******************************************************
%
% Inputs: 
%   x: a structure with fields 'lh' and 'rh', each a vector of values to
%   plot on left and right hemisphere, respectively.
%
%   atlas: atlas e.g. 'Schaefer200Yeo17', 'Glasser360', 'Yeo7', ...
%
%   Name-Value Pair Arguments (optional):
%
%   WhichOrient: 'top-bottom' or 'left-right'
%
%   MatchColorbar: logical [default: true].
% 
%   The other Name-Value Pair Arguments are as in hb_surfplot.m
%
%
% Examples:
%
%
% Dependencies:
%
%
% Hamid Behjat

d_func = fileparts(mfilename('fullpath'));
d = fullfile(d_func,'external','freesurfer','fsaverage_hb');
if exist(d, 'dir')
    d_fsavghb = d;
else
    d_fsavghb = [];
end
d = fullfile(d_func,'external','freesurfer','matlab');
if exist(d, 'dir')
    d_fsmatlab = d;
else
    d_fsmatlab = [];
end
d = fullfile(d_func,'external','plotSurfaceROIBoundary_modified');
if exist(d, 'dir')
    d_psroib = d;
else
    d_psroib = [];
end

clmap = struct;
clmap.lh = turbo(256);
clmap.rh = turbo(256);

d = inputParser;
addParameter(d,'FSavgDir', d_fsavghb);
addParameter(d,'FSMatlabDir', d_fsmatlab);
addParameter(d,'PSROIBDir', d_psroib);
addParameter(d,'BoundaryType', 'midpoint');
addParameter(d,'BoundaryAtlas', []);
addParameter(d,'BoundaryEdgeColor', [0 0 0]);
addParameter(d,'BoundaryEdgeLineWidth', 1);
addParameter(d,'WhichSurf', 'inflated'); 
addParameter(d,'SkipAssertionChecks', false);
addParameter(d,'PlotLabel', 'none');
addParameter(d,'Colormap', clmap);
addParameter(d,'OneColorbar', false);
addParameter(d,'ColorbarTickLabels', []); % applicable if ColorbarTicks is a vector
addParameter(d,'ColorbarTicks', []); % a vector, 'none', or [] ([]: default ticks)
addParameter(d,'FigureHandle', []);
addParameter(d,'DoublePlot', []);
addParameter(d,'DataRange', []);
addParameter(d,'CamLight', []);
addParameter(d,'CamView', []); 
addParameter(d,'WhichOrient', 'left-right');
addParameter(d,'MatchColorbar', true);
addParameter(d,'ShowLeftRightHemisphereTitle', true);
addParameter(d,'Boundary', []);
parse(d,varargin{:});
opts = d.Results;

if isempty(opts.FigureHandle)
    opts.FigureHandle = figure;
    set(opts.FigureHandle,'Position',[461 462 650 600]);
end

if isempty(opts.Boundary)
    opts.Boundary.lh = [];
    opts.Boundary.rh = [];
end

switch opts.WhichOrient
    case 'top-bottom'
        SBP1 = 2;
        SBP2 = 1;
    case 'left-right'
        SBP1 = 1;
        SBP2 = 2;
end

if not(isstruct(opts.Colormap))
    if ismatrix(opts.Colormap)
        d = opts.Colormap;
        opts.Colormap = struct;
        opts.Colormap.lh = d;
        opts.Colormap.rh = d;
    end
end

subplot(SBP1,SBP2,1)
axis off;
if opts.OneColorbar
    assert(ischar(opts.PlotLabel),'PlotLabel format not recognized.');
    plbl = opts.PlotLabel;
else
    if ischar(opts.PlotLabel)
        if strcmp(opts.PlotLabel, 'none')
            plbl = 'none';
        else
            if opts.ShowLeftRightHemisphereTitle
                plbl = opts.PlotLabel;
            else
                plbl = [opts.PlotLabel, ' [lh]'];
            end
        end
    else
        assert(isstruct(opts.PlotLabel),'PlotLabel format not recognized.');
        assert(isfield(opts.PlotLabel,'lh'),'lh field missing.');
        plbl = opts.PlotLabel.lh;
    end
end
[h{1}, d] = hb_surfplot(x.lh,'lh',atlas,...
    'PlotLabel',plbl,...
    'Colormap',opts.Colormap.lh,...
    'WhichSurf',opts.WhichSurf,...
    'BoundaryType',opts.BoundaryType,...
    'BoundaryAtlas',opts.BoundaryAtlas,...
    'FigureHandle',opts.FigureHandle,...
    'DoublePlot',opts.WhichOrient,...
    'OneColorbar',opts.OneColorbar,...
    'ColorbarTicks',opts.ColorbarTicks,...
    'ColorbarTickLabels',opts.ColorbarTickLabels,...
    'CamLight', opts.CamLight,...
    'CamView', opts.CamView,...
    'FSavgDir', opts.FSavgDir,...
    'FSMatlabDir', opts.FSMatlabDir,...
    'DataRange', opts.DataRange,...
    'Boundary',opts.Boundary.lh);
if isempty(opts.Boundary.lh)
    opts.Boundary.lh = d;
end

subplot(SBP1,SBP2,2)
axis off;
if opts.OneColorbar
    assert(ischar(opts.PlotLabel),'PlotLabel format not recognized.');
    plbl = opts.PlotLabel;
else
    if ischar(opts.PlotLabel)
        if strcmp(opts.PlotLabel, 'none')
            plbl = 'none';
        else
            if opts.ShowLeftRightHemisphereTitle
                plbl = opts.PlotLabel;
            else
                plbl = [opts.PlotLabel, ' [rh]'];
            end
        end
    else
        assert(isstruct(opts.PlotLabel),'PlotLabel format not recognized.');
        assert(isfield(opts.PlotLabel, 'rh'),'rh field missing.');
        plbl = opts.PlotLabel.rh;
    end
end
[h{2}, d] = hb_surfplot(x.rh,'rh',atlas,...
    'PlotLabel',plbl,...
    'Colormap',opts.Colormap.rh,...
    'WhichSurf',opts.WhichSurf,...
    'BoundaryType',opts.BoundaryType,...
    'BoundaryAtlas',opts.BoundaryAtlas,...
    'FigureHandle',opts.FigureHandle,...
    'DoublePlot',opts.WhichOrient,...
    'OneColorbar',opts.OneColorbar,...
    'ColorbarTicks',opts.ColorbarTicks,...
    'ColorbarTickLabels',opts.ColorbarTickLabels,...
    'CamLight', opts.CamLight,...
    'CamView', opts.CamView,...
    'FSavgDir', opts.FSavgDir,...
    'FSMatlabDir', opts.FSMatlabDir,...
    'DataRange',opts.DataRange,...
    'Boundary',opts.Boundary.rh);
if isempty(opts.Boundary.rh)
    opts.Boundary.rh = d;
end

if any(x.lh<0) || any(x.rh<0)
    if opts.OneColorbar
        h{1}.c.Visible = 'off';
    end
else
    if opts.MatchColorbar
        cl = zeros(1,2);
        cl(1) = min(h{1}.c.Limits(1), h{2}.c.Limits(1));
        cl(2) = max(h{1}.c.Limits(2), h{2}.c.Limits(2));
        if opts.OneColorbar
            h{1}.c.Visible = 'off';
        else
            set(h{1}.c, 'Limits', cl);
        end
        set(h{2}.c, 'Limits', cl);
    end
end
end
