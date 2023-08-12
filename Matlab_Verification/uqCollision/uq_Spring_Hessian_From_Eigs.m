function [final] = uq_Spring_Hessian_From_Eigs(X)
%UQ_SPRING_HESSIAN_FROM_EIGS gets uq energy Hessian from eigenanalysis
e0p = X(:,1);
e1p = X(:,2);
e2Perp = X(:,3);
u = norm(e2Perp);

z3 = zeros(3,1);

gu = [z3;z3;e2Perp/u];

vec0 = [z3;z3;e0p];
vec0 = vec0/norm(vec0);

vec1 = [z3;z3;e1p];
vec1 = vec1/norm(vec1);


%using mu = ep = 1 for simplicity
dPsidu = 2*(u - 1);
d2Psidu2 = 2;

eig12 = 1/u*dPsidu;

final = d2Psidu2*gu*gu' + eig12*(vec0*vec0' + vec1*vec1');
end