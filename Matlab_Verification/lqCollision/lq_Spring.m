function [psi] = lq_Spring(X)
%LQ_SPRING spring energy
%   getting some values from X
e0perp = X(1:3);
e1perp = X(4:6);
e2perp = X(7:9);
len = dot(e2perp,cross(e0perp,e1perp))/norm(cross(e0perp,e1perp));
% mu frozen at 1, epsilon frozen at 1 just for simplicity
psi = (len-1)^2;
end

