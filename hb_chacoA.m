function hb_chacoA(A,txtFileName,maxDegree,sS)
%  Creates .txt file for input to Chaco.
%  Each row i shows the adjacent vertices of vertex i
%  The .txt file will be save in the given address: txtFile_graph.
%
% Inputs:
%   A: graph adjacency matrix; can be in two forms:
%   outputTxt: full address of output .txt to be saved. 
%   maxDegree: (opt) max nodal degree; e.g. 26, 98, 124. 
%
% Hamid Behjat

if ~exist('maxDegree','var')
    maxDegree = [];
end
if ~exist('sS','var') || isempty(sS)
    sS = false;
end

Nv = size(A,2);
Ne = nnz(A)/2;

edgeWs = unique(A(find(A))); %#ok<FNDSB>

if numel(edgeWs)==1 && edgeWs(1)==1
    A_type = 'binary';
elseif numel(edgeWs)>1
    A_type = 'weighthed';
end

% Create .txt file describing the graph to Chaco
txtFile = fopen(txtFileName,'wt');

% number of veritces & edges (necessary for Chaco)
if strcmp(A_type,'binary')
    % Nv Ne
    fprintf(txtFile,'%d\t%d',Nv,Ne);
else
    % Nv Ne 1 [see Chaco manual p.20-21]
    fprintf(txtFile,'%d\t%d\t1',Nv,Ne);
end
fprintf(txtFile,'\n');

% Write edges for each vertex in one row.
for iV = 1:Nv
    
    AiV = A(:,iV);
    if isempty(maxDegree)
        indF=find(AiV);
    else
        indF=find(AiV,maxDegree);
    end
    
    if strcmp(A_type,'binary')
        d = indF'; % neighbor1 neighbor 2 ...
        fprintf(txtFile,'%d\t',d);
    else
        d = [indF'; full(AiV(indF)')]; 
        d = d(:)'; % neighbor1 edge-weight1 neighbor 2 edge-weight2 ...
        fprintf(txtFile,'%d\t',d);
    end
    fprintf(txtFile,'\n');
    
    if sS && any([iV==1,rem(iV,500)==0,iV==Nv])
        prg(iV,Nv);
    end
    
end
fclose(txtFile);
end

function prg(n,N)
l = length(num2str(N));
if n==1
    fprintf('\n Preparing graph in Chaco format: ');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end
