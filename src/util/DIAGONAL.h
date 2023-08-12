#ifndef DIAGONAL_H
#define DIAGONAL_H

#include "SETTINGS.h"

class DIAGONAL
{
public:
  DIAGONAL(const SPARSE_MATRIX& A);
  virtual VECTOR apply(const VECTOR& rhs) const;

protected:
  VECTOR _invDiag;
};

#endif
