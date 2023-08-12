#ifndef PCG_H
#define PCG_H

#include "SETTINGS.h"
#include "DIAGONAL.h"

class PCG
{
public:
  PCG(const SPARSE_MATRIX& A, const DIAGONAL& diagonal);
  ~PCG();

  VECTOR solve(const VECTOR& rhs);
  VECTOR solveAMGCL(const VECTOR& rhs);
  VECTOR solveEigenStyle(const VECTOR& rhs);
  VECTOR solveCR(const VECTOR& rhs);
  VECTOR solvePCR(const VECTOR& rhs);
  const int& iterations() { return _iterations; };
  const REAL& error()     { return _error; };

protected:
  const SPARSE_MATRIX& _A;
  const DIAGONAL& _diagonal;
    
  double _tolerance;
  int _maxIterations;
  VECTOR _x0; 

  int _iterations;
  REAL _error;
};

#endif
