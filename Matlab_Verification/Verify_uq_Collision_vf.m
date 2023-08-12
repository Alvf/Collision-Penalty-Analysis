addpath('./util');
addpath('./Projection')
addpath('./uqCollision');
addpath('./X_Conversions');
materialName = sprintf('uq Collisions (vertex-face)');

fprintf('============================================\n');
fprintf('Running numerical tests for %s\n', materialName);
fprintf('============================================\n');
fprintf('Numerically verifying expression for PK1:\n');
Verify_PK1_C(@uq_Spring, @uq_Spring_PK1_Full);
fprintf('Numerically verifying expression for dPsi/dx (with simplified frame derivative):\n');
Verify_Position_Derivative(@uq_Spring, @uq_Spring_PK1_Full, @C_vf, @dC_vfdq_simpl);
fprintf('Numerically verifying expression for Hessian:\n');
Verify_Hessian_C(@uq_Spring_PK1_Full, @uq_Spring_Hessian);
fprintf('Numerically verifying Hessian built from eigenanalysis\n');
Verify_Hessian_C(@uq_Spring_PK1_Full, @uq_Spring_Hessian_From_Eigs);
fprintf('Numerically verifying positional Hessian for uq (using simplified frame derivative and del correction):\n');
Verify_Position_Hessian_delCorrection(@uq_Spring_PK1_Full,@uq_Spring_Hessian,@C_vf,@dC_vfdq_simpl,@del_dC_vfdq);

