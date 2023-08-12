#ifndef STRAND_DIAGONAL_H
#define STRAND_DIAGONAL_H

#include "SETTINGS.h"
#include <vector>
#include "DIAGONAL.h"

class STRAND_DIAGONAL : public DIAGONAL
{
public:
  STRAND_DIAGONAL(const SPARSE_MATRIX& A, const VECTORI& strandEnds);
  virtual VECTOR apply(const VECTOR& rhs) const override;

  virtual ~STRAND_DIAGONAL();

protected:
  std::vector<SPARSE_MATRIX> _blocks;

  // even if the global matrix is SPD, it doesn't mean the blocks will be.
  // because off-diagonal entries that push the diagonal back to zero
  // have been cut out 
  //Eigen::SimplicialLDLT<SPARSE_MATRIX>* _blockInverses;
  Eigen::SparseLU<SPARSE_MATRIX>* _blockInverses;

  const VECTORI& _strandEnds;
  const int _totalStrands;
  VECTORI _strandStarts;
};

#endif
