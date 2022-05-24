function R=hb_rotmatrix(v1,v2)
% Rotation matrix mapping v1 to v2, such that v2=R*v1.
%
%
% Hamid Behjat

d = v1'*v2;
c = cross(v1,v2);
nb = norm(v1);
nd = norm(v2);
nc = norm(c);
if ~all(c==0)
    z = [0 -c(3) c(2); c(3) 0 -c(1); -c(2) c(1) 0];
    R = (eye(3)+z+z^2*(1-d)/nc^2)/nb^2;
else % collinear
    R = sign(d)*(nd/nb) ;
end
end