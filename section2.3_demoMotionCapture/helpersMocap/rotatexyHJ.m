function Vrot = rotatexyHJ(V, ang)
% rotate xy components of 3D vectors in V by angle ang 
%
% V: matrix size N x 3
% ang: angle in rad. Can be a single scalar or vector size N x 1
%
% the rotation orientation is counterclockwise in viewing from above
Vrot = V;
[theta, r] = cart2pol(Vrot(:,1), Vrot(:,2));
[Vrot1, Vrot2] = pol2cart(theta + ang, r);
Vrot(:,1:2) = [Vrot1, Vrot2];
