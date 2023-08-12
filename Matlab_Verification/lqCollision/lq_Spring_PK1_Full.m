function [pk1] = lq_Spring_PK1_Full(X)
%LQ_SPRING_PK1_Full Gets dPsi/dX without premature simplification
%(usable for finite differences)
e0 = X(1:3);
e1 = X(4:6);
e2 = X(7:9);
len = dot(e2,cross(e0,e1))/norm(cross(e0,e1));
cross01 = cross(e0,e1);
cross12 = cross(e1,e2);
cross20 = cross(e2,e0);

dlde0 = (e1*dot(e0,e1)/norm(e1)^2 - e0)*dot(e0,cross12)*norm(e1)^2/norm(cross01)^3 + cross12/norm(cross01);
dlde1 = (e0*dot(e1,e0)/norm(e0)^2 - e1)*dot(e0,cross12)*norm(e0)^2/norm(cross01)^3 + cross20/norm(cross01);
dlde2 = cross01/norm(cross01);
dldX = [dlde0',dlde1',dlde2'];
%using mu = 1 and eps = 1
dPsidlen = 2*(len-1);
pk1 = dPsidlen*dldX;
end

