function [Y, outs] = hb_graph_shuffle_v0(X,varargin)
% rewires input graph. 
% 
% Types: 
%   PreserveDegreeSequence (default); edge weights and degree distribution
%   is retained; i.e, edges are shuffled around such that the distribution
%   of binary degrees of the output graph is the same as the input graph.
%
% Inputs:
%   X: adjacency matrix of a non-directed graph with no self loops; binary
%   or weighted edges.
%
% Outputs:
%   Y: shuffled graph adjacency matrix.
%   D_in: vertex degrees input graph.
%   D_out: vertex degrees output graph.
%
% HB

d = inputParser;
addParameter(d,'Type', 'PreserveDegreeSequence');
addParameter(d,'DebugMode', false);
addParameter(d,'Verbose', false);
addParameter(d,'MaxSwap', []);
parse(d,varargin{:});
opts = d.Results;

if isempty(opts.MaxSwap)
    opts.MaxSwap = floor(0.5*0.5*nnz(X)); % half the total number of edges
end

N = size(X,1);

assert(issymmetric(X), 'input is not symmetric');
assert(isequal(diag(X),zeros(N,1)), 'non-zero diagonal (self loops)');

switch opts.Type
    case 'PreserveDegreeSequence'
        [Y, Y_only_rewire,d] = subfcn_randomize_graph_hb(X, opts.MaxSwap, opts.DebugMode);
        % Y : rewire + weight-swap on edges not rewired
        % Y_only_rewire: only rewire (no weight-swap)
        outs = struct;
        outs.Y_only_rewire = Y_only_rewire;
        outs.n_edgepairs_maxswap = opts.MaxSwap;
        outs.n_edgepairs_rewire = d.n_edgepairs_rewire;
        outs.n_edgepairs_weightswap = d.n_edgepairs_weightswap;
end
outs.MaxSwap = opts.MaxSwap;

%-Checks.
assert(issymmetric(Y), 'fishy: output not symmetric'); % symmetric?
assert(isequal(diag(Y),zeros(N,1)), 'fishy: output has self loops'); % no self-loops?
end

%==========================================================================
function e = pickedge(E,i,j,a,b)
Ii = and(i~=a, i~=b);
Ij = and(j~=a, j~=b);
Iij = and(Ii, Ij);
I = and(logical(E(:)), Iij);
d = E(I);
e = d(randi(length(d)));
% ~30% faster than:
% d = find(I);
% e = E(d(randi(length(d))));
end

%==========================================================================
function [A, A_only_rewire, out] = subfcn_randomize_graph_hb(A,MaxSwap,DebugMode)
% ----
% Based on: SUBFCN_RANDOMIZE_GRAPH in fcn_randomize_str.m 
% ----
%
% SUBFCN_RANDOMIZE_GRAPH_HB swaps edges/weights such that input degree
% sequence is preserved.
%
%   A = FCN_RANDOMIZE_GRAPH_HB(A,MaxSwap) takes adjacency matrix A and
%       performs complete rewiring, returning a randomized matrix A.
%
%   Inputs:       A,      adjacency matrix
%           MaxSwap,      maximum number of edge-pairs to swap
%         DebugMode,      binary; if true, figure shown + stats printed. 
%   Outputs:      A,      randomized matrix (entails edge/weight swaps)
%     A_only_rewire,      randomized matrix (entails only edge swaps)
%
%   Richard Betzel, Indiana University, 2013; H Behjat, 2025.
%
%   Notes:
%   1. Based on script written by Olaf Sporns (randmio_und.m).
%   2. Graph may become disconnected as a result of rewiring. Always
%      important to check.
%   3. A can be weighted, though the weighted degree sequence will not be
%      preserved.
%   4. A must be undirected.
%
%   5. [H Behjat]: SUBFCN_RANDOMIZE_GRAPH.m extended to incorporate a
%   graeter extent of randomisation while preserving the desgree sequence.
%   This is done by not only rewiring (i.e. swapping edges) but also
%   swaping edge weights. Swaping edge weights is a practical solution for
%   non-sparse graphs as otherwise there is very little options available
%   to rewire, i.e., to generate a new edge between two nodes that does not
%   exist in the input graph. In such a scenario (non-sparse, very dense
%   graphs), no matter what MaxSwap value is used, the rewred graph does
%   not differ much from the input graph, whereas resorting to "weight
%   swaping" doe the trick, which can be intuitively still considered as
%   "resiring" since existing wires are swaped.

A0 = A;
MaxEdgePairFindAttemps = 1;
D_in = sum(logical(A),1);
[i,j] = find(triu(A,1));
m = length(i);
E = 1:m; % edge indices 
Ernd1 = randperm(m,m);
nrewire = 0;

%-Step 1: swap edges, i.e. rewire. 
%--------------------------------------------------------------------------
for ie1=1:m
    e1 = Ernd1(ie1);
    if E(e1)==0
        continue; % edge already rewired
    end
    a = i(e1);
    b = j(e1);
    for ie2=1:MaxEdgePairFindAttemps
        e2 = pickedge(E,i,j,a,b);
        assert(E(e2)~=0, 'fishy');
        c = i(e2);
        d = j(e2);
        assert(all(a~=[c,d]) && all(b~=[c,d]));
        d1 = not(or(A(a,d), A(c,b)));
        d2 = not(or(A(a,c), A(d,b)));
        if d1 && d2
            if rand>0.5
                swaptype = 'ad-cb';
            else
                swaptype = 'ac-db';
            end
        elseif d1
            swaptype = 'ad-cb';
        elseif d2
            swaptype = 'ac-db';
        else
            swaptype = 'none';
        end
        switch swaptype
            case {'ad-cb', 'ac-db'}
                % Check 2 passed; no edge bw a node of e1 & e2.
                EdgePairFound = true;
            case 'none'
                % One of verices of e1 is already connected to a vertex of
                % e2, thus, rewiring as below is not possible. Note that
                % this is a major issue for close-to-complete/complete
                % graphs wherein many edges can not be rewired to a
                % candidate node that has an empty synapase to wire a new
                % edge to.
                EdgePairFound = false;
                if nnz(E)==2
                    break;
                else
                    continue;
                end
        end
        switch swaptype
            case {'ad-cb', 'ac-db'}
                switch swaptype
                    case 'ad-cb'
                        ab_new = A(a,d);
                        cd_new = A(c,b);
                        A(a,d) = A(a,b);
                        A(d,a) = A(b,a);
                        A(c,b) = A(c,d);
                        A(b,c) = A(d,c);
                    case 'ac-db'
                        ab_new = A(a,c);
                        cd_new = A(d,b);
                        A(a,c) = A(a,b);
                        A(c,a) = A(b,a);
                        A(d,b) = A(c,d);
                        A(b,d) = A(d,c);
                end
                assert(ab_new==0);
                assert(cd_new==0);
                A(a,b) = ab_new; % 0;
                A(b,a) = ab_new; % 0;
                A(c,d) = cd_new; % 0;
                A(d,c) = cd_new; % 0;
                assert(nnz(A)==nnz(A0), 'fishy');
        end
        break;
    end
    if EdgePairFound
        nrewire = nrewire + 2;
        E(e1) = 0; % mark edge as rewired
        E(e2) = 0; % mark edge as treated
        if (nrewire/2)==MaxSwap
            break;
        end
    else
        if nnz(E)==2
            break;
        end
    end
end

A_only_rewire = A;

%-Step 2: swap weights of exiting edges.
%--------------------------------------------------------------------------
% For edges that could not be rewired (which are a lot when working with
% not-sparse graph), their weights are swaped to ensure the structure of
% the input graph is destroyed, such that only the degree sequence
% preserved. Note that swapping weights can also be seen as rewiring just
% that no "new" edge is created but the existing ones are swapped in pairs.
nweightswap = MaxSwap-(nrewire/2);
if nweightswap<=0
    fprintf('\n.MaxSwap alreadty reached via rewiring; edge-weight swapping skipped.');
else
    Eskip = logical(E);
    Iskip = find(Eskip); % edges not rewired
    N0 = length(Iskip);
    assert((nweightswap*2)<=N0);
    Iskip = Iskip(randperm(N0,nweightswap*2));
    e1 = Iskip(1:nweightswap);
    e2 = Iskip(nweightswap+1:(2*nweightswap));
    a = i(e1);
    b = j(e1);
    c = i(e2);
    d = j(e2);
    sz = size(A);
    Iab = sub2ind(sz,a,b);
    Iba = sub2ind(sz,b,a);
    Icd = sub2ind(sz,c,d);
    Idc = sub2ind(sz,d,c);
    w1 = A(Iab);
    w2 = A(Icd);
    assert(isequal(A(Iba), w1));
    assert(isequal(A(Idc), w2));
    A(Iab) = w2;
    A(Iba) = w2;
    A(Icd) = w1;
    A(Idc) = w1;
end

out = struct;
out.n_edgepairs_rewire     = nrewire/2;   % number of edge-pairs rewired
out.n_edgepairs_weightswap = nweightswap; % number of edges for which only their weights were swaped

if DebugMode
    nskip = m-nrewire;
    fprintf('\n.Number of requested edge-pair swaps      : %05d (%05.1f%% of total)', MaxSwap, ((2*MaxSwap)/m)*100);
    fprintf('\n..swaps done via edge-pair rewiring       : %05d (%05.1f%% of total)', nrewire/2, (nrewire/m)*100);
    fprintf('\n..swaps done via edge-weight-pair swapping: %05d (%05.1f%% of total)\n', nweightswap, ((2*nweightswap)/m)*100);
    D_out = sum(logical(A),1);
    D_out_only_rewire = sum(logical(A_only_rewire),1);
    assert(isequal(D_out, D_in), 'fishy');
    assert(isequal(D_out_only_rewire, D_in), 'fishy');
    hf = figure;
    set(hf, 'Position', [10 100 800 500]);
    subplot(2,3,[1 2 3]);
    stem(D_in, 'Color', [0 0 1], 'MarkerSize',12);
    hold on;
    stem(D_out_only_rewire, 'fill', 'Color', [0 1 0], 'MarkerSize', 6);
    stem(D_out, 'fill', 'Color', [0 0 0], 'MarkerSize', 3);
    legend({'input', 'only rewire', 'rewire+swap'}, 'FontSize', 10);
    xlabel('nodes');
    ylabel('degree');
    title('degree sequence');
    set(gca, 'FontSize',12);
    subplot(2,3,4);
    imagesc(A0);
    axis image;
    xlabel('nodes');
    ylabel('nodes');
    title('input graph');
    set(gca, 'FontSize',12);
    subplot(2,3,5);
    imagesc(A_only_rewire);
    axis image;
    xlabel('nodes');
    ylabel('nodes');
    title('rewired graph (only rewire)');
    colormap pink;
    set(gca, 'FontSize',12);
    subplot(2,3,6);
    imagesc(A);
    axis image;
    xlabel('nodes');
    ylabel('nodes');
    title('rewired graph (rewire+swap)');
    colormap pink;
    set(gca, 'FontSize',12);
end
end
