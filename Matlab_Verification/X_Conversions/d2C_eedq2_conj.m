function [final] = d2C_eedq2_conj(pk1,numstate)
%D2C_EEDQ2_CONJ Gets pk1:d2Cee/dx2 at numstate ASSUMING pk1=0;0;nonnzero stack

state = sym('x%d',[12 1]);
assume(state, 'real');

v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1 - v0;
e1 = v3 - v2;
e2 = v2 - v0;

betSym = (e1'*e1 * e2'*e0 - e2'*e1 * e0'*e1)/(e0'*e0 * e1'*e1 - e0'*e1 * e0'*e1);
gamSym = (e0'*e0 * e2'*e1 - e2'*e0 * e0'*e1)/(e0'*e0 * e1'*e1 - e0'*e1 * e0'*e1);

dbetdx = sym('dbetdx',[12 1]);
dgamdx = sym('dgamdx', [12 1]);

for i = 1:12
    dbetdx(i) = diff(betSym, state(i));
    dgamdx(i) = diff(gamSym, state(i));
end

dbetdv0 = dbetdx(1:3);
dgamdv0 = dgamdx(1:3);

dbetdv1 = dbetdx(4:6);
dgamdv1 = dgamdx(4:6);

dbetdv2 = dbetdx(7:9);
dgamdv2 = dgamdx(7:9);

dbetdv3 = dbetdx(10:12);
dgamdv3 = dgamdx(10:12);

final = zeros(12,12);
nontriv = pk1(7:9);

%the delta portion
for i = 1:12
    b0 = diff(-e0*dbetdv0' - e1*dgamdv0', state(i));
    b1 = diff(-e0*dbetdv1' - e1*dgamdv1', state(i));
    b2 = diff(-e0*dbetdv2' - e1*dgamdv2', state(i));
    b3 = diff(-e0*dbetdv3' - e1*dgamdv3', state(i));

    b0num = double(subs(b0,state,numstate));
    b1num = double(subs(b1,state,numstate));
    b2num = double(subs(b2,state,numstate));
    b3num = double(subs(b3,state,numstate));

    final(i,:) = [nontriv'*b0num, nontriv'*b1num, nontriv'*b2num, nontriv'*b3num];
end

dbetdxnum = double(subs(dbetdx,state,numstate));
dgamdxnum = double(subs(dgamdx,state,numstate));

%the simple terms
for i = 1:12
    final(i,:) = final(i,:) + [dbetdxnum(i)*nontriv', -dbetdxnum(i)*nontriv', dgamdxnum(i)*nontriv', -dgamdxnum(i)*nontriv'];
end

end

