function opts = hb_surfplot_bothhemi(x,atlas,varargin)
% HB_SURFPLOT_BOTHHEMI same as HB_SURFPLOT.m but showing projection of data
% on both left and right hemipheres in one figure.
%
% *******************************************************
% ***** "preliminary" help; not accurate & thorough *****
% *******************************************************
%
% Inputs: 
%   x: either a vector of length equal to the number of regions in the
%   atlas, or a structure with fields 'lh' and 'rh', each a vector of
%   values to plot on left and right hemisphere, respectively.
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
addParameter(d,'FigurePosition', []);
addParameter(d,'FigureHandle', []);
addParameter(d,'FigureName', []);
addParameter(d,'DoublePlot', []);
addParameter(d,'DataRange', []);
addParameter(d,'CamLight', 'none');
addParameter(d,'CamView', []); 
addParameter(d,'WhichOrient', 'left-right');
addParameter(d,'MatchColorbar', false); % not very ueful as the colorbar just gets streched for one or both ends rather than the colors in the surface renders changing  
addParameter(d,'ShowLeftRightHemisphereTitle', true);
addParameter(d,'ShowDataRangeInFigureName', false);
addParameter(d,'ColorbarTickLabels', []); % applicable if ColorbarTicks is a vector
addParameter(d,'ColorbarTicks', []); % a 1x2 vector, 'none', or [] ([]: default ticks)
addParameter(d,'ColorbarDecimalDigits', 2);
addParameter(d,'ColorbarNumberOfTicks', 5);
addParameter(d,'Boundary', []);
parse(d,varargin{:});
opts = d.Results;

if opts.OneColorbar
    opts.MatchColorbar = true;
end

switch class(x)
    case 'struct'

    case {'double', 'logical'}
        N = floor(length(x)/2);
        d = x;
        x = struct;
        x.lh = d(1:N);
        x.rh = d(N+1:end);
end

if isempty(opts.FigureHandle)
    opts.FigureHandle = figure;
    switch opts.WhichOrient
        case 'one-row'
            if isempty(opts.FigurePosition)
                d = [50 400 1000 240];
            else
                d = opts.FigurePosition;
            end
        case {'top-bottom', 'left-right'}
            if isempty(opts.FigurePosition)
                d = [461 462 650 600];
            else
                d = opts.FigurePosition;
            end
    end
    set(opts.FigureHandle, 'Position', d);
end

if not(isempty(opts.FigureName))
    set(opts.FigureHandle, 'Name', opts.FigureName);
end

if isempty(opts.Boundary)
    opts.Boundary.lh = [];
    opts.Boundary.rh = [];
end

switch opts.WhichOrient
    case 'one-row'
        SBP1 = 1;
        SBP2 = 2;
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
    'BoundaryEdgeColor', opts.BoundaryEdgeColor,...
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
    'BoundaryEdgeColor', opts.BoundaryEdgeColor,...
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

if isempty(atlas)
    ParcelPLot = false;
else
    ParcelPLot = true;
end

if isempty(opts.ColorbarTicks)
    if ParcelPLot
        chk1 = any(x.lh<0) || any(x.rh<0);
        chk2 = any(x.lh>0) || any(x.rh>0);
    else
        chk1 = any(x.lh.data<0) || any(x.rh.data<0);
        chk2 = any(x.lh.data>0) || any(x.rh.data>0);
    end

    if chk1 && chk2
        % signed values
        d1 = max(abs(h{1}.c.Limits));
        d2 = max(abs(h{2}.c.Limits));
        cl = max(d1,d2)*[-1 1];
    else
        % all positive or negative values
        cl = zeros(1,2);
        cl(1) = min(h{1}.c.Limits(1), h{2}.c.Limits(1));
        cl(2) = max(h{1}.c.Limits(2), h{2}.c.Limits(2));
    end

    if ~isempty(opts.ColorbarDecimalDigits)
        cl = round(cl, opts.ColorbarDecimalDigits);
    end

    if isempty(opts.ColorbarDecimalDigits)
        chk1 = h{1}.c.Limits(1) == h{2}.c.Limits(1);
        chk2 = h{1}.c.Limits(2) == h{2}.c.Limits(2);
    else
        d = opts.ColorbarDecimalDigits;
        chk1 = round(h{1}.c.Limits(1),d) == round(h{2}.c.Limits(1),d);
        chk2 = round(h{1}.c.Limits(2),d) == round(h{2}.c.Limits(2),d);
    end
    OneColorbarOk = and(chk1, chk2);

    if opts.MatchColorbar
        if not(OneColorbarOk)
            warning('Not useful to MatchColorbar as lh & rh ranges differ');
            % This will just squeeze/strech the colorbar and one or both ends
            % for lh &/or rh, without affecting the colors in the surface
            % rendering. This will make it ambigious as to what the minimum or
            % maximum value is for each plot given the saturation that takes
            % place.
        end
        set(h{1}.c, 'Limits', cl);
        set(h{2}.c, 'Limits', cl);
    end

    if isempty(opts.ColorbarNumberOfTicks)
        % not applicable if two extremes of colorbar are being set
        if opts.ShowDataRangeInFigureName
            if opts.MatchColorbar
                d = sprintf('%s [%0.02f %0.02f]', ...
                    get(opts.FigureHandle, 'Name'), ...
                    cl(1), ...
                    cl(2));
                set(opts.FigureHandle, 'Name', d);
            else
                d = sprintf('%s lh:[%0.02f %0.02f] rh:[%0.02f %0.02f]', ...
                    get(opts.FigureHandle, 'Name'), ...
                    h{1}.c.Limits(1), ...
                    h{1}.c.Limits(2), ...
                    h{2}.c.Limits(1), ...
                    h{2}.c.Limits(2));
                set(opts.FigureHandle, 'Name', d);
            end
        end
    else
        % colorbar value at two extrems set
        % + ColorbarNumberOfTicks-2 more ticks in between
        d = opts.ColorbarNumberOfTicks;
        d1_a0 = h{1}.c.Limits(1);
        d1_b0 = h{1}.c.Limits(2);
        d2_a0 = h{2}.c.Limits(1);
        d2_b0 = h{2}.c.Limits(2);
        if isempty(opts.ColorbarDecimalDigits)
            d1_a = d1_a0;
            d1_b = d1_b0;
            d2_a = d2_a0;
            d2_b = d2_b0;
            h{1}.c.Ticks = linspace(d1_a, d1_b, d);
            h{2}.c.Ticks = linspace(d2_a, d2_b, d);
        else
            n = opts.ColorbarDecimalDigits;
            d1_a = fixval(d1_a0, n, 'lower');
            d1_b = fixval(d1_b0, n, 'upper');
            d2_a = fixval(d2_a0, n, 'lower');
            d2_b = fixval(d2_b0, n, 'upper');
            d1_ab = unique(round(linspace(d1_a, d1_b, d), n));
            d2_ab = unique(round(linspace(d2_a, d2_b, d), n));
            h{1}.c.Ticks = d1_ab;
            h{2}.c.Ticks = d2_ab;
        end
    end
else
    OneColorbarOk = true;
end

if opts.OneColorbar
    assert(OneColorbarOk, ...
        'Misleading to use one colorbar as the lh & rh range differ');
    % Despite MatchColorbar being done, the change in the colorbard
    % (streched or squized from one or both ends) does not translate to the
    % surface plot, and the resulting colorbars for lh and rh show
    % different values for the same color shown in lh and rh. Therefore, it
    % is not correct to use one colorbar.
    %
    % If this is not clear, you may place a break point on the error line
    % above to see how the two colorbars look like, and why it is not
    % correct to remove one of them.
    h{1}.c.Visible = 'off';
end
end

%==========================================================================
function y = fixval(x,n,WhichLim)
switch n
    case 1
        ds = 0.1;
    case 2
        ds = 0.01;
    case 3
        ds = 0.001;
end
y = round(x,n);
switch WhichLim
    case 'lower'
        if y<x
            y = y + ds;
        end
    case 'upper'
        if y>x
            y = y - ds;
        end
end
end