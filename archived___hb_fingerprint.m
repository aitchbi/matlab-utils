function [C, M] = hb_fingerprint(A,B)
% HB_FINGERPRINT return matrix and/or related measures. 
% 
% [C, M] = hb_get_identifiability(A,B)
% A: feature vectors test;    one column per subejct.
% B: feature vectors re-test; one column per subejct.
% C: identifiability matrix (IM)
% M: identifiability measures
%
% [~, measures] = hb_get_identifiability(A)
% A: identifiability matrix (IM)
% M: identifiability measures
% 
% Hamid Behjat 

N = size(A,2);

if nargin==1
    C  = A;
elseif nargin==2
    C = zeros(N, N);
    for i1=1:N
        for i2=1:N
            d = corrcoef(A(:,i1), B(:,i2));
            C(i1,i2) = d(1,2);
        end
        prgs(i1,N);
    end
end

M = struct;
M.Iself   = mean(C(find(eye(N))));
M.Iothers = mean(C(find(not(eye(N)))));
M.Idiff   = M.Iself - M.Iothers;
M.Iself_individual = diag(C)';
M.Idiff_individual = zeros(1,N);
M.accuracy         = false(1,N);
M.accuracy_rank    = zeros(1,N);

for k=1:N
    d = setdiff(1:N,k);
    d = mean(C(k,d));
    M.Idiff_individual(k) = M.Iself_individual(k) - d;
    [~, d] = sort(C(k,:), 'descend');
    M.accuracy_rank(k) = find(d==k); % subject rank
    M.accuracy(k)      = M.accuracy_rank(k)==1; % subject detected?
end

end

%==========================================================================
function prgs(n,N)
l = numel(num2str(N));
if n==1
    fprintf('\nBuilding IM..');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end
