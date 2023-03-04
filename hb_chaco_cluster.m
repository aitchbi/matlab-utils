function [C,d1,d2] = hb_chaco_cluster(A,indices,ref,path_chaco,saveDir,N,pV,maxDegree)
% HB_chaco_cluster performs spectral clustering to parcel a given graph to a
% desired number of clusteres using CHACO.
%
% Inputs:
%   A: adjacency matrix; variable or file address.
%   ref: reference volume; graph mask.
%   indices: graph indices in ref
%   path_chaco: path to Chaco executable.
%   saveDir: directort to save files in.
%   N: number of clusters.
%   pV: (opt) volume of parcels; N is determined. 
%   maxDegree: (opt) max nodal degree in graph; for speed up.
%
% Outputs: 
%   C: clustered volume; file address.
%
% Dependencies: 
%   hb_chaco_graph.m
%   hb_chaco2Vol.m
%   SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12)
%
% Note: SPM12 used only for transforming Chaco output clusters to volume.
%
% Hamid Behjat 

if ~exist('pV','var')
    pV = [];
end
if ~exist('maxDegree','var')
    maxDegree = [];
end

% Number of parcels.
assert(xor(isempty(N),isempty(pV)),'Either N or pV should be [].')
if isempty(N)
    N = round(length(indices)/pV);
end

if ischar(A)
    d = load(A);
    A = d.A;
end
assert(length(indices)==size(A,1));

% Chaco settings.
chacoSets = [2,2,1,N,1];

% File names.
f_graph = fullfile(saveDir,'A_chaco.graph.txt');                       % A [.txt]
f_chaco = fullfile(saveDir,sprintf('chaco_inputs_%04dregions.txt',N)); % chaco inputs [.txt]
f_clustT = fullfile(saveDir,sprintf('clust_%04dregions.txt',N));       % chaco output [.txt]
f_clustN = strrep(f_clustT,'.txt','.nii');                             % chaco output transformed to .nii

% Chaco inputs.
fprintf('\n.Preparing Chaco inputs.. ')
hb_chaco_graph(A,f_graph,maxDegree); % graph in Chaco format.
d = fopen(f_chaco,'wt'); 
fprintf(d,[f_graph,'\n']);
fprintf(d,[f_clustT,'\n']);
fprintf(d,[num2str(chacoSets(1)),'\n']);
fprintf(d,[num2str(chacoSets(2)),'\n']);
fprintf(d,[num2str(chacoSets(3)),'\n']);
fprintf(d,[num2str(chacoSets(4)),'\n']);
fprintf(d,[num2str(chacoSets(5)),'\n']);
fprintf(d,'n');
fclose(d);
fprintf('done.')

% Run Chaco.
fprintf('\n.Running Chaco..\n')
[d1,d2] = system([path_chaco,filesep,'chaco < ',f_chaco]);
fprintf('.done.')

% Tranform Chaco output to volume. 
fprintf('\n.Transforming Chaco output to .nii.. ')
d = hb_chaco2vol(ref,f_clustT,N,f_clustN,indices,chacoSets,[],pV);
sanityCheck(d,N,A);
fprintf('done.')

C = f_clustN;
end

function sanityCheck(I,N,A)
assert(numel(I)==N); % number of clusteres should be as expected?
for iN = 1:N % clusters connected based on A?
    dm = I{iN};
    [d,~,d2] = unique(conncomp(graph(A(dm,dm))));
    if length(d)~=1
        warning('Seems chaco hasn''t clustered as expected.');
        sz = zeros(size(d));
        for iD=1:length(d)
            sz(iD) = nnz(d2==iD);
        end
        fprintf('\n[cluter-%03d] size of conn componnets: %d',sz);
    end
end
end
