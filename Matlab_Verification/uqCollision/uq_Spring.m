function [psi] = uq_Spring(X)
%uq_SPRING Gets sqrtvf spring energy from C_matrix
e2Perp = X(:,3);
u = sqrt(e2Perp'*e2Perp);

%treating mu as 1 and ep as 1
psi = (u - 1)^2;
end

