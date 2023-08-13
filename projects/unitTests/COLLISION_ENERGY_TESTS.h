//////////////////////////////////////////////////////////////////////////////
// do a convergence test on a generic x-based collision energy
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionEnergyGradient(const ENERGY_12D* energy, 
                                            const vector<VECTOR3>& v)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING gradient for " << energy->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = energy->psi(v);
  VECTOR12 gradient = energy->gradient(v);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiffGradient;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> vNew = v;
        vNew[j][i] += eps;
        vector<VECTOR3> vDotNew = vNew;

        // get the new psi
        double psi = energy->psi(vNew);

        // store the finite difference
        finiteDiffGradient[entry] = (psi - psi0) / eps;
      }

    VECTOR12 diff = gradient - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / gradient.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " gradient: " << endl << gradient << endl;
      cout << " finite diff: " << endl << finiteDiffGradient << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test for a collision frame 
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionFrameGradient(const C_PLANES* frameFunc, 
                                            const VECTOR12& x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING gradient for " << frameFunc->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX9x12 dcdx = frameFunc->dCdx_s(x) + frameFunc->del_dCdx_s(x);
  VECTOR9 c0 = frameFunc->getC(x);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX9x12 finiteDiffGradient;

    // for each of the degrees of the freedom
    for (int i = 0; i < 12; i++)
    {
      VECTOR12 xNew = x;
      xNew[i] += eps;

      // get the new c
      VECTOR9 c = frameFunc->getC(xNew);

      // store the finite difference
      finiteDiffGradient.block<9,1>(0,i) = (c - c0) / eps;
    }

    MATRIX9x12 diff = dcdx - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / dcdx.norm())) / (9*12.0);
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " dcdx:" << endl << dcdx << endl;
      cout << " finite diff:" << endl << finiteDiffGradient << endl;
      cout << " diff:" << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on an x-based Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionEnergyHessian(const ENERGY_12D* energy, 
                                           const vector<VECTOR3>& v)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Hessian for " << energy->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  VECTOR12 gradient0 = energy->gradient(v);
  MATRIX12 H = energy->hessian(v);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiffH;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> vNew = v;
        vNew[j][i] += eps;

        // get the new psi
        VECTOR12 gradient = energy->gradient(vNew);

        // store the finite difference
        finiteDiffH.col(entry) = (gradient - gradient0) / eps;
      }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// test out collision energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Collision Energy tests", "[collision]" )
{
  // so we test against consistent examples
  reseed(123456);

  vector<VECTOR3> vertices(4);
  for (int x = 0; x < 4; x++)
    vertices[x] = randomVector3();
  VECTOR12 x = randomVector12();
  
  const REAL mu = 1000.0;
  const REAL collisionEps = 0.1;
  SECTION("Testing signed vertex-face collision energy")
  {
    ENERGY_1D normalSpring(mu, collisionEps);
    SIGNED_LEN_PLANES lengthFunc;
    C_PLANES cFunc;
    
    ENERGY_12D vertexFaceEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&vertexFaceEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&vertexFaceEnergy, vertices) == true);
  }

  SECTION("Testing cubic signed vertex-face collision energy")
  {
    CUBIC_SIGNED_SPRING_1D normalSpring(mu, collisionEps);
    SIGNED_LEN_PLANES lengthFunc;
    C_PLANES cFunc;
    
    ENERGY_12D vertexFaceEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&vertexFaceEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&vertexFaceEnergy, vertices) == true);
  }
  
  SECTION("Testing unsigned vertex-face collision energy")
  {
    UNSIGNED_SPRING_1D normalSpring(mu, collisionEps);
    UNSIGNED_LEN_PLANES lengthFunc;
    C_PLANES cFunc;
    
    ENERGY_12D vertexFaceEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&vertexFaceEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&vertexFaceEnergy, vertices) == true);
  }
  
  SECTION("Testing unsigned cubic spring vertex-face collision energy")
  {
    CUBIC_SIGNED_SPRING_1D normalSpring(mu, collisionEps);
    UNSIGNED_LEN_PLANES lengthFunc;
    C_PLANES cFunc;
    
    ENERGY_12D vertexFaceEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&vertexFaceEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&vertexFaceEnergy, vertices) == true);
  }

  SECTION("Testing signed edge-edge collision energy")
  {
    ENERGY_1D normalSpring(mu, collisionEps);
    SIGNED_LEN_PLANES lengthFunc;
    C_PLANES_EE cFunc;
    
    ENERGY_12D edgeEdgeEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&edgeEdgeEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&edgeEdgeEnergy, vertices) == true);
  }

  SECTION("Testing signed cubic spring edge-edge collision energy")
  {
    CUBIC_SIGNED_SPRING_1D normalSpring(mu, collisionEps);
    SIGNED_LEN_PLANES lengthFunc;
    C_PLANES_EE cFunc;
    
    ENERGY_12D edgeEdgeEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&edgeEdgeEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&edgeEdgeEnergy, vertices) == true);
  }

  SECTION("Testing unsigned edge-edge collision energy")
  {
    UNSIGNED_SPRING_1D normalSpring(mu, collisionEps);
    UNSIGNED_LEN_PLANES lengthFunc;
    C_PLANES_EE cFunc;
    
    ENERGY_12D edgeEdgeEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&edgeEdgeEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&edgeEdgeEnergy, vertices) == true);
  }

  SECTION("Testing unsigned cubic spring edge-edge collision energy")
  {
    CUBIC_SIGNED_SPRING_1D normalSpring(mu, collisionEps);
    UNSIGNED_LEN_PLANES lengthFunc;
    C_PLANES_EE cFunc;
    
    ENERGY_12D edgeEdgeEnergy(&normalSpring, &cFunc, &lengthFunc);
    REQUIRE(convergenceTestCollisionEnergyGradient(&edgeEdgeEnergy, vertices) == true);
    REQUIRE(convergenceTestCollisionEnergyHessian(&edgeEdgeEnergy, vertices) == true);
  }

  SECTION("Testing eeSpecial")
  {
  C_SPECIAL_EE cFunc;
  REQUIRE(convergenceTestCollisionFrameGradient(&cFunc,x) == true);
  }
}

