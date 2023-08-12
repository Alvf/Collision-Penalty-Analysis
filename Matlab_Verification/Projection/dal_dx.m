function [dadx] = dal_dx(state)
%DAL_DX gets deriv of al for collision projection frame

v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1 - v2;
e1 = v3 - v2;

e0h = e0/norm(e0);

dadx = 1/(e0'*e0) * [ zeros(3,1);
                      e1 - 2*e0h*(e0h'*e1);
                      2*e0h*(e0h'*e1) - e1 - e0;
                      e0];

end

