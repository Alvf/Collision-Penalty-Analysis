function t = Verify_Hessian_C(PK1,Hessian)
%VERIFY_HESSIAN_X Verifying simplified Hessian (expecting orthogonal frames
%   random diagonal
S = 10*[rand,0,0;0,rand,0;0,0,rand];
%   random rotation
R = quat2rotm(randrot);
C = R*S;
% track whether or not the differences are converging
previousdiff = inf;
failed = false;  
  
Hdirect = Hessian(C);
Hnumerical = zeros(9,9);
eps = 1e-2;
  for e = 1:5
    X0 = C;
    P0 = PK1(X0);
    
    i = 1;
    for y = 1:3
      for x = 1:3
        X1 = C;
        X1(x,y) = X1(x,y) + eps;
        P1 = PK1(X1);
        diff = (P1 - P0) / eps;
        
        Hnumerical(:,i) = reshape(diff, 9,1);
        i = i + 1;
      end
    end
  
    Hdiff = Hdirect - Hnumerical;
    diff = norm(Hdiff) / 81;
    
    % see if we stopped converging
    if (diff >= (previousdiff / 2))
        failed = true;
    end
    previousdiff = diff;
    
    fprintf('eps: %.10f \t diff: %.10f\n', eps, diff);
    eps = eps * 0.1;
  end
  if (failed)
      fprintf('Hessian convergence test ***FAILED***\n')
      Hdirect
      Hnumerical
      [V,D] = eig(Hnumerical)
  else
      fprintf('Hessian convergence test ***PASSED***\n')
  end  
end

