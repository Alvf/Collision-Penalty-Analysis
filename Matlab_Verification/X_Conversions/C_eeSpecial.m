function [final] = C_eeSpecial(state)
%C_EESPECIAL Gets special frame for edge-edge uq case
%   get vertices
v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1-v0;
e1 = v3-v2;
e2 = v2-v0;

% normally a and b would be taken from finding the vector of
% least distance, but since they're just treated as constant
% we can, without loss of generality, set a and b to some constant
% and trust generality (if you change a and b here make sure to also 
% change it inside dC_eeSpecialdx_simpl)

a = 0.25;
b = 0.46;

final = [e0, e1, e2 - a*e0 + b*e1];
end

