addpath('./Gen_Length');
addpath('./X_Conversions');

%Numerically verifying that A1 + A2 = -A1 in the paper (H_l case with vertex-face)

%get a random starting state
state = rand(12,1);

fprintf('Getting A1:\n');
C = C_vf(state);

delTerm = del_dC_vfdq(state);
Hl = Gen_Length_Hessian(C);
A1 = delTerm'*Hl*delTerm

fprintf('Getting A2:\n');
simpleTerm = dC_vfdq_simpl(state);
gl = reshape(Gen_Length_PK1_Full(C),9,1);

contractionTerm = d2C_vfdq2_conj(gl, state);
A2 = delTerm'*Hl*simpleTerm + simpleTerm'*Hl*delTerm + contractionTerm

fprintf('Getting A1 + A2:\n')
A1pA2 = A1 + A2

fprintf('Comparing to -A1 (should be effectively zero):\n')
A1pA1pA2 = A1 + A1pA2