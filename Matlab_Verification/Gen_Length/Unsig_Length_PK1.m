function [pk1] = Unsig_Length_PK1(C)
%UNSIG_LENGTH_PK1 Gets pk1 of u starting from given C
%   Detailed explanation goes here
e2p = C*[0;0;1];
u = norm(e2p);
pk1 = 1/u*[zeros(3,1), zeros(3,1), e2p];
end

