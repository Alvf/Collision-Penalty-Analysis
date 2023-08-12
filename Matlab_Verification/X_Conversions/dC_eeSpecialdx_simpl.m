function [simpl] = dC_eeSpecialdx_simpl(state)
%DC_EESPECIALDX_SIMPL Gets (dCeeSpecial/dx)_s.

%don't even technically need state vector

% normally a and b would be taken from finding the vector of
% least distance, but since they're just treated as constant
% we can, without loss of generality, set a and b to some constant
% and trust generality (if you change a and b here make sure to also 
% change it inside C_eeSpecial)

a = 0.25;
b = 0.46;

id = eye(3);
z3 = zeros(3);
simpl = [-id, id, z3,z3;
         z3, z3, -id, id;
         a*id - id, -a*id, -b*id+id, b*id];
end

