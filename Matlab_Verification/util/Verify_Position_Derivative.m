function t = Verify_Position_Derivative(Psi, PK1, get_C, dCdq)

%Generate a random starting state
state = rand(12,1);

%Get X from state
C = get_C(state);

%calculate the analytic position derivative
derivativeDirect = dCdq(state)' * reshape(PK1(C),9,1);

derivativeNumerical = zeros(12,1);
eps = 1e-2;

% track whether or not the differences are converging
previousdiff = inf;
failed = false;

for iter = 1:5
    state0 = state;
    psi0 = Psi(get_C(state0));
    
    for i = 1:12
        state1 = state0;
        state1(i) = state1(i) + eps;
        psi1 = Psi(get_C(state1));
        diff = (psi1 - psi0)/eps;

        derivativeNumerical(i) = diff;
    end
    
    derivativeDiff = derivativeDirect - derivativeNumerical;
    diff = norm(derivativeDiff)/12;

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
    derivativeDirect
    derivativeNumerical
else
    fprintf('Derivative convergence test ***PASSED***\n')
end

end