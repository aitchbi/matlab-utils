% [Hamid Behjat]
% Based on plotSurfaceROIBoundary.m; cleaned up & added following options:
% .option to change edgecolor from the default(black).
% .


function [p,boundary_plot,BOUNDARY] = plotSurfaceROIBoundary_hb(data,vertex_id,surface,opts) %cmap,varargin)

% This script is a wrapper for findROIboundaries and makeFaceVertexCData so
% that they smoothly work together and you don't have to spend a lot of
% your own time making them work together. The code sets up the basics of
% the patch object

% Inputs:
%
% surface = a structure with two fields: vertices (the vertices making up
% the surface) and faces (the faces of the surface)
%
% vertex_id = the roi id of each vertex
%
% data = either data for each individual roi or data for each vertex.
% If you don't want any data to be displayed for a roi or vertex, set that
% value to NaN. Note that this assumes a roi with a vertex_id has no data
% to plot. Additionally, if the vertex_ids are non sequential (e.g., like
% what you get from an .annot file) then data can take the form of a ROI*2
% matrix, where ROI is the number of regions, each row is a particular
% region with the first column being the data to plot and the second being
% the region ID (should correspond to what is in vertex_id)
%
% boundary_method = 'faces', 'midpoint', 'centroid', 'edge_vertices', or
% 'edge_faces'. 'faces' will find the faces which exist between ROIs and
% those will be coloured black to specify the boundary. 'midpoint' finds
% the edges that connect the vertices of two different ROIs and takes the
% midpoint of the edge and uses those coordinates to define the boundary.
% 'centroid' finds the faces which exist between ROIs and uses the centroid
% of those to draw the coordinates that define the boundary.
% 'edge_vertices' finds the vertices which define the boundary of the ROI
% and uses them for the coordinates. 'edge_faces' finds the edges along
% faces which make up the boundary and uses them for the coordinates
%
% cmap = an N*3 matrix specifying the RGB values making up the colormap to
% use
%
% linewidth = the width of the boundary when using 'midpoint', or
% 'centroid'.
%
% 'ColorbarLimits' = the range to apply the colormap. This will work perfectly fine
% if the range specified is larger than the data itself or if the upper
% limit is larger. However if the lower limit is larger than the smallest
% value in data, then it may get strange. If colorUnknownGrey = 1, then
% faces/vertices with a value smaller than the limit will be coloured grey,
% or potentially black if 'faces' is used. If it is set to 0 and 'faces' is
% used, those regions will be set to black while if 'centroid' or midpoint'
% are selected, the colormap will work appropriately). So if you really
% need to enforce a lower limit I would suggest threshold the data in
% advance and all should be good.
%
% Outputs:
%
% p = the patch surface object
%
% boundary_plot = a strcture containing a line object defining each
% boundary
%
% BOUNDARY = the boundary of the rois. For 'faces' this will be a logical
% where a value of 1 indicates that face is on the boundary between ROIs.
% For 'midpoint', 'centroid' or 'edges', BOUNDARY will be a cell where each
% element contains the coordinates of the points making up the boundary,
% which can be plotted as a line. Note that each boundary defines a
% continuous ROI, if a ROI is made up of multiple non-continuous parts
% (i.e., the ROI is made up of multiple unconnected sections), there will
% be a boundary for each of those parts


%d = inputParser;
%addParameter(d,'BoundaryType',[]);
%addParameter(d,'BoundaryEdgeColor',[0 0 0]);
%addParameter(d,'BoundaryEdgeLineWidth',2);
%addParameter(d,'ColorbarLimits',[]);
%addParameter(d,'Boundary',[]);
%parse(d,varargin{:});
%opts = d.Results;

% Extract the faces and vertices
vertices = surface.vertices;
faces = surface.faces;

assert(length(vertices)==length(vertex_id));

if length(data)==length(vertices)
    PlotResolution = 'vertex';
else
    PlotResolution = 'parcel';
end

assert(isvector(data));
% data, either:
% .vector with a value per ROI
% .vector with a value per vertex

data = data(:);

switch PlotResolution
    case 'parcel'
        N_labels = unique(vertex_id);
        if isequal(N_labels(:), [0:3, 5:35]')
            % label 4 (corpuscallosum) missing for Desikan-Killiany (aparc) atlas, but ok.
            N_labels = length(N_labels)+1;
        else
            N_labels = length(N_labels);
        end

        if any(vertex_id==0)
            % ignore label 0 if present, which is medial wall
            N_labels = N_labels-1;
        end
        assert(length(data)==N_labels);

    case 'vertex'
        assert(length(data)==length(vertices));
end

switch opts.BoundaryType
    case 'faces'
        ColorFaceBoundaries = true;
        FaceColor = 'flat'; % has to be 'flat' for patch to work

    case {'midpoint', 'centroid', 'edge_vertices', 'edge_faces'}

        ColorFaceBoundaries = false;
        FaceColor = 'flat';
        assert(strcmp(PlotResolution, 'parcel'));

    case 'none'
        ColorFaceBoundaries = false;
        FaceColor = 'interp';
end

%-Plot surface map.
%--------------------------------------------------------------------------
% 1. plot patch
p = patch(surface);

% 2. assign colors
FVCD = ...
    getFVCD(...
    data,...
    vertex_id, ...
    faces,...
    opts.Colormap,...
    opts.DataRange,...
    PlotResolution,...
    ColorFaceBoundaries);
set(p,...
    'FaceVertexCData',FVCD,...
    'FaceColor',FaceColor,...
    'EdgeColor','none',...
    'Clipping','off');

p.FaceLighting = 'gouraud';

material dull

if strcmp(PlotResolution, 'vertex')
    BOUNDARY = [];
    boundary_plot = [];
    return;
end

%-Find parcel boundaries.
%--------------------------------------------------------------------------
if isempty(opts.Boundary)
    switch opts.BoundaryType
        case {'midpoint', 'centroid', 'edge_vertices', 'edge_faces', 'faces'}
            BOUNDARY = findROIboundaries( ...
                vertices, ...
                faces, ...
                vertex_id, ...
                opts.BoundaryType);

        case 'none'
            BOUNDARY = [];
    end
else
    BOUNDARY = opts.Boundary;
end

%-Draw parcel boundaries.
%--------------------------------------------------------------------------
switch opts.BoundaryType
    case {'midpoint', 'centroid', 'edge_vertices', 'edge_faces'}
        hold on
        for k = 1:length(BOUNDARY)
            boundary_plot.boundary(k) = plot3(...
                BOUNDARY{k}(:,1),...
                BOUNDARY{k}(:,2),...
                BOUNDARY{k}(:,3),...
                'Color', opts.BoundaryEdgeColor,...
                'LineWidth',opts.BoundaryEdgeLineWidth,...
                'Clipping','off');
        end

    case 'none'
        boundary_plot = [];
end
end

%==========================================================================
function FVCD = getFVCD(data,vertex_id,faces,cmap,DataRange,PlotRes,CFB,varargin)
%**************************************************************************
% Based on makeFaceVertexCData.m by Stuart Oldham, Monash University, 2020.
%**************************************************************************
% This script will plot the boundaries defined by some parcellation on
% a surface projection. It can also alter the colormap so regions that do
% not have any information are coloured grey.
%
% ColorLimits: [--USE WITH CAUTION!--] the range to apply the colormap.
% This will work perfectly fine if the range specified is larger than the
% data itself or if the upper limit is larger. However if the lower limit
% is larger than the smallest value in data, then it may get strange. If
% colorUnknownGrey = 1, then faces/vertices with a value smaller than the
% limit will be coloured grey. If it is set to 0 and 'faces' is used, those
% regions will be set to black (if 'centroid' or midpoint' are selected,
% the colormap will work appropriately). So if you really need to enforce a
% lower limit I would suggest threshold the data in advance and all should
% be good.
%
% ColorFaceBoundaries: logical. If true, faces that make up boundaries of
% each parcel are coloured (color: ColorForBoundary). If true, the surface
% plot will have a value per face instead of a value per vertex.
%
% ColorForUnknown: color for unknown regions (vertex_id==0 or data==NaN).
%
% ColorForBoundary: color for boundary if drawn.
%
% Output:
% FVCD = the color value for each vertex or face (depending on how
% colorFaceBoundaries was configured).

d = inputParser;
addParameter(d,'ColorForUnknown', [.5 .5 .5]); % unknown parcels/vertices
addParameter(d,'ColorForBoundary', [0 0 0]);
parse(d,varargin{:});
opts = d.Results;

N_data = length(data);

switch PlotRes
    case 'parcel'
        N_parcel = N_data;
        ParcelID = 1:N_parcel;
        
    case 'vertex'
        % n/a
end

if isempty(DataRange)
    switch PlotRes
        case 'parcel'
            cmin = nanmin(data);
            cmax = nanmax(data);
        case 'vertex'
            cmin = nanmin(data(vertex_id>0));
            cmax = nanmax(data(vertex_id>0));
    end
else
    assert(all(not(isnan(DataRange))));
    cmin = DataRange(1);
    cmax = DataRange(2);
end

switch PlotRes
    case 'parcel'
         newval = [NaN; data(:)];
         oldval = [0; ParcelID(:)];

    case 'vertex'
        % n/a
end

if CFB
    faces_roi_ids = vertex_id(faces); % parcels each face is connected to
    face_roi_id = faces_roi_ids(:,1); % parcel id of each face
    FaceBoundary = logical(diff(faces_roi_ids,2,2)); % boundary faces
    switch PlotRes
        case 'parcel'
            fv_data = face_roi_id;
            for k = 1:(N_parcel+1)
                fv_data(face_roi_id == oldval(k)) = newval(k);
            end
            % face value: value of parcel that the face belongs to

        case 'vertex'
            fv_data = nanmean(data(faces),2);
            fv_data(face_roi_id==0) = NaN;
            % face value: mean of values of its associated vertices
    end
else
    FaceBoundary = [];
    switch PlotRes
        case 'parcel'
            fv_data = vertex_id;
            for k = 1:(N_parcel+1)
                fv_data(vertex_id == oldval(k)) = newval(k);
            end

        case 'vertex'
            fv_data = data;
            fv_data(vertex_id==0) = NaN;
    end
end

FVCD = assigncolors(fv_data, cmap, cmin, cmax, opts, CFB, FaceBoundary);
end

%=======================================================================
function FVCD = assigncolors(fv_data,cmap,cmin,cmax,opts,CFB, FB)

assert(size(cmap,2)==3);

N_data = length(fv_data);
N_cmap = size(cmap,1);
I_nan  = isnan(fv_data);

if 1 
    % new/accurate approach
    
    % Map to an index for the colormap
    d = rescale(fv_data, 1, N_cmap, 'InputMin', cmin, 'InputMax', cmax);
    color_ind = round(d);
    
    assert(isequal(isnan(color_ind), I_nan));

    color_ind(I_nan) = 1; % temporary color for NaNs (unknown regions)

else % old/messy/fishy approach

    fv_data(fv_data<cmin) = cmin;

    fv_data(fv_data>cmax) = cmax;

    fv_data(N_data+1) = cmin;
    fv_data(N_data+2) = cmax;

    % Map to an index for the colormap
    color_ind = ceil(rescale(fv_data, 1, size(cmap,1)));

    asert(isequal(isnan(color_ind), I_nan));

    % Temporarily assign NaNs (i.e., the value for unknown regions)
    % to a value so logical indexing doesn't screw up
    color_ind(isnan(color_ind)) = 1;

end

% assign color to faces
FVCD = cmap(color_ind(1:N_data),:);

% assign color to unknown regions
FVCD(I_nan,1) = opts.ColorForUnknown(1);
FVCD(I_nan,2) = opts.ColorForUnknown(2);
FVCD(I_nan,3) = opts.ColorForUnknown(3);

% assign color to boundary faces
if CFB
    FVCD(FB,1) = opts.ColorForBoundary(1);
    FVCD(FB,2) = opts.ColorForBoundary(2);
    FVCD(FB,3) = opts.ColorForBoundary(3);
end

%FVCD = FVCD'; % 3 x number of faces
end
