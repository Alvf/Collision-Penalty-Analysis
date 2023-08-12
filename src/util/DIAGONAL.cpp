#include "DIAGONAL.h"

DIAGONAL::DIAGONAL(const SPARSE_MATRIX& A)
{
  _invDiag.resize(A.rows());

  for (unsigned int x = 0; x < A.rows(); x++)
    _invDiag[x] = 1.0;

  // see which blocks exist
  for (unsigned int x = 0; x < A.outerSize(); x++)
    for (SPARSE_MATRIX::InnerIterator it(A, x); it; ++it)
    {
      //if (it.row() == it.col())
      if (it.index() == x)
        _invDiag[it.row()] = 1.0 / it.value();
    }
}

VECTOR DIAGONAL::apply(const VECTOR& rhs) const
{
  VECTOR result(rhs.size());
  
  for (unsigned int x = 0; x < result.size(); x++)
    result[x] = _invDiag[x] * rhs[x];
    
  return result;
}
