function [final] = Gen_Length_Hessian(X)
%GEN_LENGTH_HESSIAN Hessian of general length
%   Recover perps from X
e0perp = X(1:3);
e1perp = X(4:6);
e2perp = X(7:9);
len = dot(e2perp,cross(e0perp,e1perp))/norm(cross(e0perp,e1perp));
zeroBlock = zeros(3,3);

e0perpHat = e0perp/norm(e0perp);
e1perpHat = e1perp/norm(e1perp);
projectionBlock = eye(3) - e0perpHat'*e0perpHat - e1perpHat'*e1perpHat;
d2lde0de2 = e2perp'*e0perp/(norm(e2perp)^2*norm(e0perp)^2);
d2lde1de2 = e2perp'*e1perp/(norm(e2perp)^2*norm(e1perp)^2);
final = -len*[projectionBlock/norm(e0perp)^2, zeroBlock, d2lde0de2;
              zeroBlock, projectionBlock/norm(e1perp)^2, d2lde1de2;
              d2lde0de2', d2lde1de2', zeroBlock];
end

