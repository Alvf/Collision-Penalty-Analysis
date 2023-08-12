function [dgdx] = dgam_dx(numstate)
%DGAM_DX gets dgamdx from collision projection frame
v0 = numstate(1:3);
v1 = numstate(4:6);
v2 = numstate(7:9);
v3 = numstate(10:12);

e0 = v1 - v2;
e1 = v3 - v2;
e2 = v0 - v2;

gam = (e0'*e0 * e2'*e1 - e2'*e0 * e0'*e1)/(e0'*e0 * e1'*e1 - e0'*e1 * e0'*e1);

e1n2 = e1'*e1;
e0n2 = e0'*e0;
e0d1 = e0'*e1; 
e0d2 = e0'*e2;
e1d2 = e1'*e2;
e0x12 = cross(e0,e1)'*cross(e0,e1);

dgdv0 = -e0d1*e0 + e0n2*e1;
dgdv1 = 2*gam*(e0d1*e1 - e1n2*e0) + 2*e1d2*e0 - e0d2*e1 - e0d1*e2;
dgdv2 = 2*gam*(e1n2*e0 - e0d1*(e0 + e1) + e0n2*e1) - 2*e1d2*e0 + e0d2*(e0+e1) + e0d1*(e0 + e2) - e0n2*(e1+e2);
dgdv3 = 2*gam*(e0d1*e0 - e0n2*e1) - e0d2*e0 + e0n2*e2;

dgdx = [dgdv0;dgdv1;dgdv2;dgdv3]/e0x12;

end

