function [hess] = Unsig_Length_Hessian(C)
%UNSIG_LENGTH_HESSIAN Gets H_u

e2p = C*[0;0;1];
u = norm(e2p);

z33 = zeros(3,3);
i3 = eye(3);

hess = [z33, z33, z33;
        z33, z33, z33;
        z33, z33, 1/u*(i3 - e2p*e2p'/(e2p'*e2p))];

end

