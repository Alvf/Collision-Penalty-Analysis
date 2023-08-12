function[final] = Gen_Length(X)
%GEN_LENGTH Returns general signed length from state frame
%   get vectors from X
e0perp = X(1:3);
e1perp = X(4:6);
e2perp = X(7:9);
final = dot(e2perp,cross(e0perp,e1perp))/norm(cross(e0perp,e1perp));
end

