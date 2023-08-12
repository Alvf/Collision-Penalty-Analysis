function [final] = Gen_Length_Hessian_From_Eigs(X)
%GEN_LENGTH_HESSIAN_FROM_EIGS Constructs Hessian from eigenpairs 
%   Recover perps from X
e0perp = X(1:3);
e1perp = X(4:6);
e2perp = X(7:9);
len = dot(e2perp,cross(e0perp,e1perp))/norm(cross(e0perp,e1perp));

zeroVec = [0;0;0];

lam0 = -len/(2*norm(e1perp)^2)*(1 + sqrt(1 + 4 * norm(e1perp)^2/norm(e2perp)^2));
al0 = lam0;
bet0 = lam0 - len/norm(e2perp)^2;
vec0 = [zeroVec;bet0*e2perp';al0*e1perp'];
vec0 = vec0/norm(vec0);

lam1 = -len/(2*norm(e1perp)^2)*(1 - sqrt(1 + 4 * norm(e1perp)^2/norm(e2perp)^2));
al1 = lam1;
bet1 = lam1 - len/norm(e2perp)^2;
vec1 = [zeroVec;bet1*e2perp';al1*e1perp'];
vec1 = vec1/norm(vec1);

lam2 = -len/(2*norm(e0perp)^2)*(1 + sqrt(1 + 4 * norm(e0perp)^2/norm(e2perp)^2));
al2 = lam2;
bet2 = lam2 - len/norm(e2perp)^2;
vec2 = [bet2*e2perp';zeroVec;al2*e0perp'];
vec2 = vec2/norm(vec2);

lam3 = -len/(2*norm(e0perp)^2)*(1 - sqrt(1 + 4 * norm(e0perp)^2/norm(e2perp)^2));
al3 = lam3;
bet3 = lam3 - len/norm(e2perp)^2;
vec3 = [bet3*e2perp';zeroVec;al3*e0perp'];
vec3 = vec3/norm(vec3);

diag = zeros(9,9);
diag(1,1) = lam0;
diag(2,2) = lam1;
diag(3,3) = lam2;
diag(4,4) = lam3;

vecs = zeros(9,9);
vecs(:,1) = vec0;
vecs(:,2) = vec1;
vecs(:,3) = vec2;
vecs(:,4) = vec3;

final = vecs*diag*vecs';

end

