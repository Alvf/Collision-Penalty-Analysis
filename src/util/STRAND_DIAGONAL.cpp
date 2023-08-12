#include "STRAND_DIAGONAL.h"
#include "TIMER.h"
#include "MATRIX_UTIL.h"

#include <iostream>

using namespace std;

STRAND_DIAGONAL::STRAND_DIAGONAL(const SPARSE_MATRIX& A, const VECTORI& strandEnds) :
  DIAGONAL(A), _strandEnds(strandEnds), _totalStrands(strandEnds.size())
{
  TIMER functionTimer("Build STRAND_DIAGONAL preconditioner");
  // track where everything starts too
  _strandStarts = strandEnds;

  // copy the blocks in
  _blocks.resize(_totalStrands);
  int start = 0;
  for (int x = 0; x < _totalStrands; x++)
  {
    const int end = strandEnds[x] - 1;
    const int length = end - start + 1;
    _strandStarts[x] = start;
    _blocks[x] = A.block(start,start,length,length);

    /*
    SPARSE_MATRIX myBlock = A.block(start,start,length,length);
    cout << " A original: " << endl << MATRIX(A) << endl;
    cout << " Block: " << endl << MATRIX(myBlock) << endl;

    cout << " start: " << start << " length: " << length << endl;
    */

    start = end + 1;
  }
  
  // factorize
  cout << " Building preconditioner ..." << flush;
  // even if the global matrix is SPD, it doesn't mean the blocks will be.
  // because off-diagonal entries that push the diagonal back to zero
  // have been cut out 
  //_blockInverses = new Eigen::SimplicialLDLT<SPARSE_MATRIX>[_totalStrands];
  _blockInverses = new Eigen::SparseLU<SPARSE_MATRIX>[_totalStrands];
#pragma omp parallel
#pragma omp for schedule(static)
  for (int x = 0; x < _totalStrands; x++)
    _blockInverses[x].compute(_blocks[x]);

  // print some info if the inverse bombed 
  for (int x = 0; x < _totalStrands; x++)
  {
    if (_blockInverses[x].info() != Eigen::Success)
    {
      cout << " Block inverse failed on strand " << x << ". Eigenvalues are:" << endl;
      cout << HOBAK::eigenvalues(_blocks[x]) << endl;
      cout << " Info: " << _blockInverses[x].info() << endl;
      exit(0); 
    }
    assert(_blockInverses[x].info() == Eigen::Success);
  }
  cout << " done" << endl;
}

STRAND_DIAGONAL::~STRAND_DIAGONAL()
{
  delete[] _blockInverses;
}

VECTOR STRAND_DIAGONAL::apply(const VECTOR& rhs) const
{
  TIMER functionTimer(string("STRAND_DIAGONAL::") + __FUNCTION__);
  VECTOR result(rhs.size());
  result.setZero();

#pragma omp parallel
#pragma omp for schedule(static)
  for (int x = 0; x < _totalStrands; x++)
  {
    const int start = _strandStarts[x];
    const int end = _strandEnds[x];
    const int length = end - start;
    //cout << " start: " << start << " end: " << end << endl;
    assert(end <= result.size());
    const VECTOR rhsSegment = rhs.segment(start, length);
    //cout << " rhs dims: " << rhsSegment.size() << endl;

    // get the appropriate segment
    //result.segment(start, end) = _blockInverses[x].solve(rhs.segment(start, end));
    const VECTOR resultSegment = _blockInverses[x].solve(rhsSegment);
    result.segment(start, length) = resultSegment;
  }

  //cout << " Diagonal solution: " << result.transpose() << endl;
    
  return result;
}
