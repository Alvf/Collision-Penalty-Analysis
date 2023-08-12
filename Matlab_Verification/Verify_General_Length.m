addpath('./util');
addpath('./Gen_Length');

fprintf('============================================\n');
fprintf('Running numerical tests for %s\n', "Signed length (l)");
fprintf('============================================\n');
fprintf('Numerically verifying expression for PK1:\n');
Verify_PK1_C(@Gen_Length, @Gen_Length_PK1);
fprintf('Numerically verifying expression for Hessian:\n');
Verify_Hessian_C(@Gen_Length_PK1_Full, @Gen_Length_Hessian);
fprintf('Numerically verifying eigenpairs for Hessian:\n');
Verify_Hessian_C(@Gen_Length_PK1_Full, @Gen_Length_Hessian_From_Eigs);

fprintf('============================================\n');
fprintf('Running numerical tests for %s\n', "Unsigned length (u)");
fprintf('============================================\n');
fprintf('Numerically verifying expression for PK1:\n');
Verify_PK1_C(@Unsig_Length, @Unsig_Length_PK1);
fprintf('Numerically verifying expression for Hessian:\n');
Verify_Hessian_C(@Unsig_Length_PK1, @Unsig_Length_Hessian);
fprintf('Numerically verifying eigenpairs for Hessian:\n');
Verify_Hessian_C(@Unsig_Length_PK1, @Unsig_Length_Hessian_From_Eigs);