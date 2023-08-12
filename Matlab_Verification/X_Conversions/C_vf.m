function [final] = C_vf(state)
%C_VF get Cvf from 12x1 state vector

%   get vertices
v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1 - v2;
e1 = v3 - v2;
e2 = v0 - v2;

e1perp = e1 - e1'*e0/(e0'*e0)*e0;

A = [e0'*e0, e0'*e1; 
     e0'*e1, e1'*e1];
b = [e2'*e0;
     e2'*e1];
betgam = A^(-1)*b;

e2perp = e2 - betgam(1)*e0 - betgam(2)*e1;

final = [e0,e1perp,e2perp];
end

