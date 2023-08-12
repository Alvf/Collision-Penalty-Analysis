function [state, ExtraTerm, HessianDirect] = Verify_Position_Hessian_delCorrection(PK1, Hessian, get_C, dCdq, deldCdq)

%Generate a random starting state
state = rand(12,1);
%state = [1;1;1; 2;0;0; 0;0;0; 1;2;0];

%Get X from state
C = get_C(state);

%calculate the analytic position derivative
HessianDirect = dCdq(state)' * Hessian(C) * dCdq(state) - deldCdq(state)'*Hessian(C)*deldCdq(state);

HessianNumerical = zeros(12,12);
eps = 1e-2;

% track whether or not the differences are converging
previousdiff = inf;
failed = false;

for iter = 1:5
    state0 = state;
    diff_0 = dCdq(state0)'*reshape(PK1(get_C(state0)),9,1);
    
    for i = 1:12
        state1 = state0;
        state1(i) = state1(i) + eps;
        diff_1 = dCdq(state1)'*reshape(PK1(get_C(state1)),9,1);
        hess_col = (diff_1 - diff_0)/eps;

        HessianNumerical(:,i) = hess_col;
    end
    
    hessDiff = HessianDirect - HessianNumerical;
    diff = norm(hessDiff)/144;

    % see if we stopped converging
    if (diff >= (previousdiff / 2))
        failed = true;
    end
    previousdiff = diff;
    
    fprintf('eps: %.10f \t diff: %.10f\n', eps, diff)
    eps = eps * 0.1;
end

if (failed)
    fprintf('Derivative convergence test ***FAILED***\n')
    HessianDirect
    [Vanal,Danal] = eigs(HessianDirect)
    HessianNumerical
    [Vnum,Dnum] = eigs(HessianNumerical)
    ExtraTerm = HessianNumerical - HessianDirect
    [vEx,dEx] = eigs(ExtraTerm)
else
    fprintf('Derivative convergence test ***PASSED***\n')
    ExtraTerm = HessianNumerical - HessianDirect;
end

end