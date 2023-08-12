function [u] = Unsig_Length(X)
%UNSIG_LENGTH Gets unsigned length from X
e2Perp = X*[0;0;1];
u = norm(e2Perp);
end

