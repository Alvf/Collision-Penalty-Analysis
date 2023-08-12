function [final] = C_ee(state)
%C_EE Gets edge edge frame
%   get vertices
v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1-v0;
e1 = v3-v2;
e2 = v2-v0;

e0x12 = cross(e0,e1)'*cross(e0,e1);

al = e0'*e1/(e0'*e0);
bet = (e1'*e1*(e0'*e2) - (e2'*e1)*(e0'*e1))/e0x12;
gam = (e0'*e0*(e2'*e1) - (e2'*e0)*(e0'*e1))/e0x12;

final = [e0, e1 - al*e0, e2 - bet*e0 - gam*e1];
%final column is equivalent to eab projected down to span(e0,e1)

end

