function t = Verify_PK1_C(Psi,PK1)
%VERIFY_PK1_X Verifies PK1s of functions expecting orthogonal frames
%   random diagonal
S = 10*[rand,0,0;0,rand,0;0,0,rand];
%   random rotation
R = quat2rotm(randrot);
C = R*S;

PK1direct = PK1(C);
PK1numerical = zeros(3,3);
eps = 1e-2;
  
  % track whether or not the differences are converging
  previousdiff = inf;
  failed = false;
  
  for e = 1:5
    X0 = C;
    P0 = Psi(X0);
    
    for y = 1:3
      for x = 1:3
        X1 = C;
        X1(x,y) = X1(x,y) + eps;
        
        P1 = Psi(X1);
        diff = (P1 - P0) / eps;
        
        PK1numerical(x,y) = diff;
      end
    end
  
    PK1diff = PK1direct - PK1numerical;
    diff = norm(PK1diff) / 9;
    
    % see if we stopped converging
    if (diff >= (previousdiff / 2))
        failed = true;
    end
    previousdiff = diff;
    
    fprintf('eps: %.10f \t diff: %.10f\n', eps, diff)
    eps = eps * 0.1;
  end
  if (failed)
      fprintf('PK1 convergence test ***FAILED***\n')
      PK1direct
      PK1numerical
  else
      fprintf('PK1 convergence test ***PASSED***\n')
  end

end

