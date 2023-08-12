addpath('./Gen_Length');
addpath('./X_Conversions');

%Numerically verifying that A1 + A2 = -A1 in the paper (H_u case with edge-edge)

%get a random starting state
state = rand(12,1);

fprintf('Getting A1:\n');
C = C_ee(state);

delTerm = del_dC_eedx(state);
Hu = Unsig_Length_Hessian(C);
A1 = delTerm'*Hu*delTerm

fprintf('Getting A2:\n');
simpleTerm = dC_eedx_simpl(state);
% gu = gl up to a sign; normally not important since it only ever shows
% up in outer products, but it does matter here.
gu = reshape(Unsig_Length_PK1(C),9,1);

contractionTerm = d2C_eedq2_conj(gu, state);
A2 = delTerm'*Hu*simpleTerm + simpleTerm'*Hu*delTerm + contractionTerm

fprintf('Getting A1 + A2:\n')
A1pA2 = A1 + A2

fprintf('Comparing to -A1 (should be effectively zero):\n')
A1pA1pA2 = A1 + A1pA2