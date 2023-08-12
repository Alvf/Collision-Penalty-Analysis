function [pk1] = uq_Spring_PK1_Full(X)
%UQ_SPRING_PK1_FULL gets full PK1 from sqrtvf
e2Perp = X(:,3);
u = norm(e2Perp);
zeroVec = zeros(3,1);
dudC = [zeroVec,zeroVec,e2Perp/u];

%treating mu as 1 and ep as 1
dPsidu = 2*(u - 1);

pk1 = dPsidu*dudC;
end

