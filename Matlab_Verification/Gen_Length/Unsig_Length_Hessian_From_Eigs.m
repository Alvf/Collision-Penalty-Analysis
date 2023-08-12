function [hess] = Unsig_Length_Hessian_From_Eigs(C)
%UNSIG_LENGTH_HESSIAN Gets H_u

e0p = C(1:3)';
e1p = C(4:6)';
e2p = C(7:9)';
u = norm(e2p);

z3 = zeros(3,1);
vec0 = [z3;z3;e0p];
vec0 = vec0/norm(vec0);
vec1 = [z3;z3;e1p];
vec1 = vec1/norm(vec1);


hess = 1/u*(vec0*vec0' + vec1*vec1');

end

