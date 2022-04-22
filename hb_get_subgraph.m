function [Asub,Isub] = hb_get_subgraph(A,I,tol,verbose)
% HB_GET_SUBGRAPH extracts a connected subgraph of a graph.
%
% A = hb_get_subgraph(A,I)
% extracts a subgraph of the graph given by adjacency matrix A whose
% vertices are those given in vector I. If teh resulting subgraph is not
% connected, it will be made so by only keeping the largest connected
% component.
%
% A = hb_get_subgraph(A,I,tol)
% tol, a percentage of nnz(I), specifies the tolerance level in ensuring
% that the extracted subgraph is connected. The default tol is 5; to skip
% check, i.e., allowing any number of vertices be removed, set to 100; to
% disallow removing any vertices, set to 0.
%
% A = hb_get_subgraph(___,verbose)
% print stats and info; default=true.

% Inputs:
%   A: graph adjacency matrix. 
%
%   I: a logcical vector of length M = size(A,1), or a double vector of
%   length <= M specifying the vertex indices of teh requested subgraph of
%   A.
%
%   tol: a scalar value within [0,100); default=5. To skip tolerance check,
%   set to 100.
% 
% Hamid Behjat.

if ~exist('tol','var') || isempty(tol)
    tol = 5;
end

if ~exist('verbose','var')
    verbose = true;
end

N = size(A,1);

% check I
assert(isvector(I));
if islogical(I)
    assert(length(I)==N);
    LogicalI = true;
else
    assert(max(I)<=N);
    d = I;
    I = false(size(A,1),1);
    I(d) = true;
    LogicalI = false;
end

% extract subgraph
Asub = A(I,I);

% check connectivity
comps = conncomp(graph(Asub));
N_comps = max(comps);
if N_comps==1
    if verbose
        fprintf('\n..subgraph conencted; no vertices removed.');
    end
    Isub = I;
    return;
end

% make Asub connected
comps_size = histcounts(comps,1:N_comps+1);
[~,d] = sort(comps_size,'descend');
Isub = ismember(comps,d(1));
d = length(Isub)-nnz(Isub);
assert(d<tol*size(Asub,1),...
    'fishy; notable number of disconnected vertices.');
if verbose
    if d==1
        s = 'vertex';
    else
        s = 'vertices';
    end
    fprintf('\n..subgraph made conencted; %d %s removed.',d,s);
end
Asub = Asub(Isub,Isub);

d = false(size(I));
d(I) = Isub;
Isub = d;

if ~LogicalI
    Isub = find(Isub);
end
end