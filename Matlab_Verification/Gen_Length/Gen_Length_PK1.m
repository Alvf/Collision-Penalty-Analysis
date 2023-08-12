function [final] = Gen_Length_PK1(X)
%GEN_LENGTH_PK1 get unflattened gl from state frame
%   Recover perps from X
e0perp = X(1:3);
e1perp = X(4:6);
e2perp = X(7:9);
len = dot(e2perp,cross(e0perp,e1perp))/norm(cross(e0perp,e1perp));
zeroVec = [0;0;0];
final = [zeroVec, zeroVec, e2perp'/norm(e2perp)^2*len];
end

