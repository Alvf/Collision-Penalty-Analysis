/*
This file is part of HOBAK.

HOBAK is free software: you can redistribute it and/or modify it under the terms of 
the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

HOBAK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with HOBAK. 
If not, see <https://www.gnu.org/licenses/>.
*/

#include "BLOCK_SPARSE_MATRIX3.h"
#include "util/TIMER.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3::BLOCK_SPARSE_MATRIX3(const int& totalRows,
                                           const int& totalCols) :
  _totalRows(totalRows), _totalCols(totalCols),
  _totalBlockRows(totalRows / 3), _totalBlockCols(totalCols / 3)
{
  _blocks.resize(_totalBlockRows);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3::BLOCK_SPARSE_MATRIX3(const SPARSE_MATRIX& A) :
  _totalRows(A.rows()), _totalCols(A.cols()),
  _totalBlockRows(A.rows() / 3), _totalBlockCols(A.cols() / 3)
{
  TIMER convertTimer("Convert to block");
  _blocks.resize(_totalBlockRows);

#if 0
  const int totalThreads = omp_get_max_threads();
  std::vector<std::unordered_map<int, bool> > blockExists(totalThreads);

  // see which blocks exist
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < A.outerSize(); x++)
  {
    const int threadID = omp_get_thread_num();

    for (SPARSE_MATRIX::InnerIterator it(A, x); it; ++it)
    {
      const int blockRow = it.row() / 3;
      const int blockCol = it.col() / 3;

      const int blockIndex = blockRow + blockCol * _totalBlockRows;
      blockExists[threadID][blockIndex] = true;
    }
  }

  // copy them into current map
  for (int x = 0; x < totalThreads; x++)
    for (auto iter = blockExists[x].begin(); iter != blockExists[x].end(); iter++)
    {
      const int blockIndex = iter->first;
      const int blockRow = blockIndex % _totalBlockRows;
      const int blockCol = blockIndex / _totalBlockRows;

      const MATRIX3 block = A.block(3 * blockRow, 3 * blockCol, 3,3);

      //_blocks[blockIndex] = block;
      _blocks[blockRow][blockCol] = block;
      _totalNonZeroBlocks++;
    }
#else
  // see which blocks exist
  for (unsigned int x = 0; x < A.outerSize(); x++)
    for (SPARSE_MATRIX::InnerIterator it(A, x); it; ++it)
    {
      const int blockRow = it.row() / 3;
      const int blockCol = it.col() / 3;

      const int blockIndex = blockRow + blockCol * _totalBlockRows;
      _blockExists[blockIndex] = true;
    }

  // copy them into current map
  for (auto iter = _blockExists.begin(); iter != _blockExists.end(); iter++)
  {
    const int blockIndex = iter->first;
    const int blockRow = blockIndex % _totalBlockRows;
    const int blockCol = blockIndex / _totalBlockRows;

    const MATRIX3 block = A.block(3 * blockRow, 3 * blockCol, 3,3);

    //_blocks[blockIndex] = block;
    _blocks[blockRow][blockCol] = block;
    _totalNonZeroBlocks++;
  }
#endif
}
 
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX3::setZero()
{
  TIMER functionTimer(__FUNCTION__);
  for (int x = 0; x < _totalBlockRows; x++)
    for (auto iter = _blocks[x].begin(); iter != _blocks[x].end(); iter++)
      iter->second = MATRIX3::Zero();
}

///////////////////////////////////////////////////////////////////////
// take the transpose
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3 BLOCK_SPARSE_MATRIX3::transpose() const
{
  const BLOCK_SPARSE_MATRIX3& current = *this;
  BLOCK_SPARSE_MATRIX3 A = current;

  for (int x = 0; x < _totalBlockRows; x++)
    for (auto iter = _blocks[x].begin(); iter != _blocks[x].end(); iter++)
    {
      //const int blockIndex = iter->first;
      //const int blockRow = blockIndex % _totalBlockRows;
      //const int blockCol = blockIndex / _totalBlockRows;
      const int blockRow = x;
      const int blockCol = iter->first;
      A(blockRow, blockCol) = iter->second.transpose();
    }

  return A;
}

///////////////////////////////////////////////////////////////////////
// access a block matrix
///////////////////////////////////////////////////////////////////////
MATRIX3& BLOCK_SPARSE_MATRIX3::operator()(const int blockRow, const int blockCol)
{
  assert(blockRow >= 0);
  assert(blockCol >= 0);
  assert(blockRow < _totalBlockRows);
  assert(blockCol < _totalBlockCols);

  // see if we created it already
  //if (_blocks.find(index) != _blocks.end())
  //  return _blocks[index];
  if (_blocks[blockRow].find(blockCol) != _blocks[blockRow].end())
    return _blocks[blockRow][blockCol];

  // if we didn't create it, initialize to zero and return it
  const int index = blockRow + blockCol * _totalBlockRows;
  _blockExists[index] = true;
  _blocks[blockRow][blockCol] = MATRIX3::Zero();
  return _blocks[blockRow][blockCol];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
const bool BLOCK_SPARSE_MATRIX3::blockExists(const int blockRow, const int blockCol) const
{
  const int blockIndex = blockRow + blockRow * _totalBlockRows;

  return !(_blockExists.find(blockIndex) == _blockExists.end());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3 operator+(const BLOCK_SPARSE_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B)
{
  BLOCK_SPARSE_MATRIX3 C = A;

  assert(A.rows() == B.rows());
  assert(A.cols() == B.cols());

  for (int x = 0; x < C._totalBlockRows; x++)
    C(x,x) = C(x,x) + B(x,x);

  return C;
}

/*
///////////////////////////////////////////////////////////////////////
// UNTESTED
///////////////////////////////////////////////////////////////////////
VECTOR operator*(const BLOCK_SPARSE_MATRIX3& A, const VECTOR& v)
{
  VECTOR y(v.size());
  y.setZero();

  for (int x = 0; x < A._totalBlockRows; x++)
    for (auto iter = A._blocks[x].begin(); iter != A._blocks[x].end(); iter++)
    {
      const int blockCol = iter->first;
      const MATRIX3& blockA = iter->second;
      const VECTOR3& blockV = v.block(3 * blockCol, 0, 3, 1);
      y.block(3 * blockCol,0,3,1) = blockA * blockV;
    }

  return y;
}
*/

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3 operator*(const BLOCK_SPARSE_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B)
{
  BLOCK_SPARSE_MATRIX3 C(A.rows(), A.cols());

  for (int x = 0; x < A._totalBlockRows; x++)
    for (auto iter = A._blocks[x].begin(); iter != A._blocks[x].end(); iter++)
    {
      //const int blockIndex = iter->first;
      //const int blockRow = blockIndex % A._totalBlockRows;
      //const int blockCol = blockIndex / A._totalBlockRows;
      const int blockRow = x;
      const int blockCol = iter->first;
      const MATRIX3& blockA = iter->second;
      const MATRIX3& blockB = B(blockCol, blockCol);
      C(blockRow, blockCol) = blockA * blockB;
    }

  return C;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3 operator*(const BLOCK_DIAGONAL_MATRIX3& A, const BLOCK_SPARSE_MATRIX3& B)
{
  BLOCK_SPARSE_MATRIX3 C(A.rows(), A.cols());
  for (int x = 0; x < B._totalBlockRows; x++)
    for (auto iter = B._blocks[x].begin(); iter != B._blocks[x].end(); iter++)
    {
      //const int blockIndex = iter->first;
      //const int blockRow = blockIndex % B._totalBlockRows;
      //const int blockCol = blockIndex / B._totalBlockRows;
      const int blockRow = x;
      const int blockCol = iter->first;
      const MATRIX3& blockA = A(blockRow, blockRow);
      const MATRIX3& blockB = iter->second;
      C(blockRow, blockCol) = blockA * blockB;
    }
  return C;
}

///////////////////////////////////////////////////////////////////////
// scale on both the left and right by matrix A
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3 BLOCK_SPARSE_MATRIX3::scaleLeftRightAdd(const BLOCK_DIAGONAL_MATRIX3& A,
                                                             const BLOCK_DIAGONAL_MATRIX3& B)
{
  BLOCK_SPARSE_MATRIX3 C(A.rows(), A.cols());
  for (int x = 0; x < _totalBlockRows; x++)
    for (auto iter = _blocks[x].begin(); iter != _blocks[x].end(); iter++)
    {
      //const int blockIndex = iter->first;
      //const int blockRow = blockIndex % _totalBlockRows;
      //const int blockCol = blockIndex / _totalBlockRows;
      const int blockRow = x;
      const int blockCol = iter->first;

      // do a symmetric multiply
      if (blockRow > blockCol) continue;

      const MATRIX3& blockRowA = A(blockRow, blockRow);
      const MATRIX3& blockRowB = B(blockRow, blockRow);

      const MATRIX3& block = iter->second;
      if (blockRow == blockCol)
      {
        const MATRIX3 product = blockRowA * block * blockRowA;
        C(blockRow, blockCol) = product + blockRowB;
        continue;
      }

      const MATRIX3& blockColA = A(blockCol, blockCol);
      const MATRIX3 product = blockRowA * block * blockColA;
      C(blockRow, blockCol) = product;
      C(blockCol, blockRow) = product.transpose();
    }
  return C;
}

///////////////////////////////////////////////////////////////////////
// scale on both the left and right by matrix A
///////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX3 BLOCK_SPARSE_MATRIX3::scaleLeftRight(const BLOCK_DIAGONAL_MATRIX3& A)
{
  BLOCK_SPARSE_MATRIX3 C(A.rows(), A.cols());
  for (int x = 0; x < _totalBlockRows; x++)
    for (auto iter = _blocks[x].begin(); iter != _blocks[x].end(); iter++)
    {
      //const int blockIndex = iter->first;
      //const int blockRow = blockIndex % _totalBlockRows;
      //const int blockCol = blockIndex / _totalBlockRows;
      const int blockRow = x;
      const int blockCol = iter->first;

      // do a symmetric multiply
      if (blockRow > blockCol) continue;

      const MATRIX3& blockRowA = A(blockRow, blockRow);
      const MATRIX3& blockB = iter->second;
      if (blockRow == blockCol)
      {
        const MATRIX3 product = blockRowA * blockB * blockRowA;
        C(blockRow, blockCol) = product;
        continue;
      }

      const MATRIX3& blockColA = A(blockCol, blockCol);
      const MATRIX3 product = blockRowA * blockB * blockColA;
      C(blockRow, blockCol) = product;
      C(blockCol, blockRow) = product.transpose();
    }
  return C;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX BLOCK_SPARSE_MATRIX3::sparseView() const
{
  TIMER functionTimer(__FUNCTION__);
  typedef Eigen::Triplet<REAL> TRIPLET;
  std::vector<TRIPLET> triplets;
  for (int i = 0; i < _totalBlockRows; i++)
    for (auto iter = _blocks[i].begin(); iter != _blocks[i].end(); iter++)
    {
      //const int blockIndex = iter->first;
      //const int blockRow = blockIndex % _totalBlockRows;
      //const int blockCol = blockIndex / _totalBlockRows;
      const int blockRow = i;
      const int blockCol = iter->first;
      const MATRIX3& block = iter->second;

      for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++)
        {
          const REAL& entry = block(x, y);

          // don't bother to pass near-zeros
          //if (fabs(entry) < 1e-7) // too severe, only gets to 42
          //if (fabs(entry) < 1e-8) // gets up to 269 
          //if (fabs(entry) < 1e-10) // gets up to 299
          if (fabs(entry) < 1e-9) // passes up to timestep 326
            continue;

          const int row = 3 * blockRow + x;
          const int col = 3 * blockCol + y;
          const TRIPLET triplet(row, col, entry);

          assert(3 * blockRow + x < _totalRows);
          assert(3 * blockCol + y < _totalCols);
          
          assert(3 * blockRow + x >= 0);
          assert(3 * blockCol + y >= 0);

          triplets.push_back(triplet);
        }
    }

  SPARSE_MATRIX C(_totalRows, _totalCols);
  C.setFromTriplets(triplets.begin(), triplets.end());
  return C;
}
