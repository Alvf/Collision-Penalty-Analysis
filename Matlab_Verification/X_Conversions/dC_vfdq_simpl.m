function [final] = dC_vfdq_simpl(state)
%DC_VFDQ_SIMPL gets (dCvf/dx)_s from the state
v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1 - v2;
e1 = v3 - v2;
e2 = v0 - v2;

alph = e1'*e0/(e0'*e0);
A = [e0'*e0, e0'*e1; 
     e0'*e1, e1'*e1];
b = [e2'*e0;
     e2'*e1];
betgam = A^(-1)*b;

Id = eye(3);
zeros3 = zeros(3);

final = [zeros3, Id, -Id, zeros3;
         zeros3, -alph*Id, (alph - 1)*Id, Id;
         Id, -betgam(1)*Id, (betgam(1) + betgam(2) - 1)*Id, -betgam(2)*Id];
end

