function [final] = uq_Spring_Hessian(X)
%UQ_SPRING_HESSIAN gets sqrtvf Hessian
e2Perp = X(:,3);
u = norm(e2Perp);
gu = [0;0;0;0;0;0;e2Perp/u];
z3 = zeros(3,3);
Hu = [z3,z3,z3;
      z3,z3,z3;
      z3,z3,1/u*(eye(3) - e2Perp*e2Perp'/u^2) ];

%using mu = ep = 1 for simplicity
dPsidu = 2*(u - 1);
d2Psidu2 = 2;

final = d2Psidu2*gu*gu' + dPsidu*Hu;
end

