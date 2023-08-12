addpath('./util');
addpath('./lqCollision');
addpath('./X_Conversions');
materialName = sprintf('lq Collisions (vertex-face)');

fprintf('============================================\n');
fprintf('Running numerical tests for %s\n', materialName);
fprintf('============================================\n');
fprintf('Numerically verifying expression for PK1:\n');
Verify_PK1_C(@lq_Spring, @lq_Spring_PK1_Full);
fprintf('Numerically verifying positional derivative for collision (using simplified frame derivative):\n');
Verify_Position_Derivative(@lq_Spring, @lq_Spring_PK1_Full, @C_vf, @dC_vfdq_simpl);
fprintf('Numerically verifying expression for Hessian:\n');
Verify_Hessian_C(@lq_Spring_PK1_Full, @lq_Spring_Hessian);
fprintf('Numerically verifying eigenpairs for Hessian:\n');
Verify_Hessian_C(@lq_Spring_PK1_Full, @lq_Spring_Hessian_From_Eigs);
fprintf('Numerically verifying positional hessian for collision (using only simplified frame derivative):\n');
Verify_Position_Hessian(@lq_Spring_PK1_Full,@lq_Spring_Hessian,@C_vf,@dC_vfdq_simpl);
