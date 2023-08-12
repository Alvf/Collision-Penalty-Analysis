function [final] = lq_Spring_Hessian(X)
%LQ_SPRING_HESSIAN lq spring energy's C-based Hessian
%   first deriv of l
e0perp = X(1:3);
e1perp = X(4:6);
e2perp = X(7:9);
len = dot(e2perp,cross(e0perp,e1perp))/norm(cross(e0perp,e1perp));
zeroVec = [0;0;0];
dldX = [zeroVec, zeroVec, e2perp'/norm(e2perp)^2*len];
gl = reshape(dldX,9,1);

%   hessian of l
zeroBlock = zeros(3,3);
e0perpHat = e0perp/norm(e0perp);
e1perpHat = e1perp/norm(e1perp);
projectionBlock = eye(3) - e0perpHat'*e0perpHat - e1perpHat'*e1perpHat;
d2lde0de2 = e2perp'*e0perp/(norm(e2perp)^2*norm(e0perp)^2);
d2lde1de2 = e2perp'*e1perp/(norm(e2perp)^2*norm(e1perp)^2);
Hl = -len*[projectionBlock/norm(e0perp)^2, zeroBlock, d2lde0de2;
              zeroBlock, projectionBlock/norm(e1perp)^2, d2lde1de2;
              d2lde0de2', d2lde1de2', zeroBlock];

%using mu = 1 and eps = 1
dPsidlen = 2*(len-1);
d2Psidlen2 = 2;

final = d2Psidlen2*(gl*gl') + dPsidlen*Hl;
end

