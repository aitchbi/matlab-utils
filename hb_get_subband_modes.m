function V = hb_get_subband_modes(U, e, g)
% HB_GET_SUBBAND_MODES buils band-specific modes.
%
% 
% h behjat

N = size(U,1);
J = length(g);
V = zeros(N,J);
for j=1:J
    w = g{j}(e(:));
    V(:,j) = U*w;
end
end