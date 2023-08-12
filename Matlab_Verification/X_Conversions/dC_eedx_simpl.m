function [simpl] = dC_eedx_simpl(state)
%DC_EEDX_SIMPL Gets (dCee/dx)_s.
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

id = eye(3);
z3 = zeros(3);
simpl = [-id, id, z3,z3;
         al*id, -al*id, -id, id;
         bet*id - id, -bet*id, gam*id+id, -gam*id];
end

