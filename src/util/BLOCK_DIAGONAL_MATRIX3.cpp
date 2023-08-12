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

#include "BLOCK_DIAGONAL_MATRIX3.h"
#include "util/TIMER.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX3::BLOCK_DIAGONAL_MATRIX3() :
  _totalRows(-1), _totalCols(-1), _totalBlockRows(-1), _totalBlockCols(-1)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX3::BLOCK_DIAGONAL_MATRIX3(const int& totalRows,
                                               const int& totalCols) :
  _totalRows(totalRows), _totalCols(totalCols),
  _totalBlockRows(totalRows / 3), _totalBlockCols(totalCols / 3)
{
  // matrix size a factor of 3, right?
  assert(totalRows % 3 == 0);
  assert(totalCols % 3 == 0);

  // only supporting square matrices for now
  assert(_totalRows == _totalCols);

  _blocks.resize(_totalBlockRows);
  for (int x = 0; x < _totalBlockRows; x++)
    _blocks[x] = MATRIX3::Zero();
}
  
///////////////////////////////////////////////////////////////////////
// access a block matrix
//
// if you try to access something off-diagonal, this will bomb
///////////////////////////////////////////////////////////////////////
MATRIX3& BLOCK_DIAGONAL_MATRIX3::operator()(const int blockRow, const int blockCol)
{
  assert(blockRow >= 0);
  assert(blockCol >= 0);
  assert(blockRow < _totalBlockRows);
  assert(blockCol < _totalBlockCols);

  // you're indexing a diagonal block, right?
  assert(blockRow == blockCol);

  return _blocks[blockRow];
}

///////////////////////////////////////////////////////////////////////
// access a block matrix
//
// if you try to access something off-diagonal, this will bomb
///////////////////////////////////////////////////////////////////////
const MATRIX3& BLOCK_DIAGONAL_MATRIX3::operator()(const int blockRow, const int blockCol) const
{
  assert(blockRow >= 0);
  assert(blockCol >= 0);
  assert(blockRow < _totalBlockRows);
  assert(blockCol < _totalBlockCols);

  // you're indexing a diagonal block, right?
  assert(blockRow == blockCol);

  return _blocks[blockRow];
}

///////////////////////////////////////////////////////////////////////
// return a transposed version
///////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX3 BLOCK_DIAGONAL_MATRIX3::transpose() const
{
  BLOCK_DIAGONAL_MATRIX3 A(_totalRows, _totalCols);

  for (int x = 0; x < _totalBlockRows; x++)
    A._blocks[x] = _blocks[x].transpose();

  return A;
}

///////////////////////////////////////////////////////////////////////
// access a scalar entry
//
// if you try to access something off-diagonal, this will bomb
///////////////////////////////////////////////////////////////////////
REAL& BLOCK_DIAGONAL_MATRIX3::entry(const int row, const int col)
{
  const int blockRow = row / 3;

  // you're indexing a diagonal block, right?
  assert(blockRow == (col / 3));
  
  const int subRow = row % 3;
  const int subCol = col % 3;

  MATRIX3& block = _blocks[blockRow];
  return block(subRow, subCol);
}

///////////////////////////////////////////////////////////////////////
// access a scalar entry
//
// if you try to access something off-diagonal, this will bomb
///////////////////////////////////////////////////////////////////////
const REAL BLOCK_DIAGONAL_MATRIX3::entry(const int row, const int col) const
{
  const int blockRow = row / 3;

  // you're indexing a diagonal block, right?
  assert(blockRow == (col / 3));
  assert(blockRow >= 0);
  assert(blockRow < _totalBlockRows);
  
  const int subRow = row % 3;
  const int subCol = col % 3;

  assert(subRow < 3);
  assert(subCol < 3);

  const MATRIX3& block = _blocks[blockRow];
  return block(subRow, subCol);
}

///////////////////////////////////////////////////////////////////////
// is there even an entry to retrieve?
///////////////////////////////////////////////////////////////////////
const bool BLOCK_DIAGONAL_MATRIX3::exists(const int row, const int col) const
{
  const int blockRow = row / 3;
  const int blockCol = col / 3;
  return blockRow == blockCol;
}

///////////////////////////////////////////////////////////////////////
// subtract two of the same type
///////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX3 operator-(const BLOCK_DIAGONAL_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B)
{
  assert(A._totalRows == B._totalRows);
  assert(A._totalCols == B._totalCols);

  BLOCK_DIAGONAL_MATRIX3 C(A._totalRows, A._totalCols);
  for (int x = 0; x < A._totalBlockRows; x++)
    C._blocks[x] = A._blocks[x] - B._blocks[x];

  return C;
}

///////////////////////////////////////////////////////////////////////
// add against a sparse matrix
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator+(const SPARSE_MATRIX& A, const BLOCK_DIAGONAL_MATRIX3& B)
{
  typedef Eigen::Triplet<REAL> TRIPLET;
  std::vector<TRIPLET> triplets;

  for (int x = 0; x < B._totalBlockRows; x++)
  {
    const MATRIX3& block = B._blocks[x];
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        const int row = 3 * x + i;
        const int col = 3 * x + j;
        
        TRIPLET triplet(row, col, block(i,j));
        triplets.push_back(triplet);
      }
  }

  SPARSE_MATRIX sparseB(B._totalRows, B._totalCols);
  sparseB.setFromTriplets(triplets.begin(), triplets.end());

  return A + sparseB;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR operator*(const BLOCK_DIAGONAL_MATRIX3& A, const VECTOR& v)
{
  VECTOR y(v.size());
  y.setZero();

  using namespace std;

  for (int x = 0; x < A._totalBlockRows; x++)
  {
    const MATRIX3& blockA = A._blocks[x];
    const VECTOR3& blockV = v.segment(3 * x, 3);
    y.segment(3 * x, 3) = blockA * blockV;
  }

  return y;
}

/*
///////////////////////////////////////////////////////////////////////
// multiply against a sparse matrix
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(const SPARSE_MATRIX& A, const BLOCK_DIAGONAL_MATRIX3& B)
{
  assert(A.rows() == B._totalRows);
  assert(A.cols() == B._totalCols);

  SPARSE_MATRIX C = B.transpose() * A.transpose();
  return C.transpose();
}
*/

//////////////////////////////////////////////////////////////////////
// block matrix - full matrix multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(const BLOCK_DIAGONAL_MATRIX3& A, const SPARSE_MATRIX& B)
{
  TIMER functionTimer(__FUNCTION__);
  assert(A._totalRows == B.rows());
  assert(A._totalCols == B.cols());
  const int& totalBlockRows = A._totalBlockRows;

  //std::map<int, bool> blockExistsB;
  std::vector<bool> blockExistsB(totalBlockRows * totalBlockRows);
  for (int x = 0; x < totalBlockRows * totalBlockRows; x++)
    blockExistsB[x] = false;
  
  // see which blocks exist
  for (unsigned int x = 0; x < B.outerSize(); x++)
    for (SPARSE_MATRIX::InnerIterator it(B, x); it; ++it)
    {
      const int blockRow = it.row() / 3;
      const int blockCol = it.col() / 3;

      const int blockIndex = blockRow + blockCol * totalBlockRows;
      blockExistsB[blockIndex] = true;
    }

  typedef Eigen::Triplet<REAL> TRIPLET;
  std::vector<TRIPLET> triplets;

  TIMER multiply("Multiply");
  // go through all the blocks
  for (int i = 0; i < totalBlockRows; i++)
  {
    const MATRIX3& blockA = A._blocks[i];
    for (int j = 0; j < totalBlockRows; j++)
      for (int k = 0; k < totalBlockRows; k++)
      {
        if (i != k) continue;

        int blockIndex = k + j * totalBlockRows;
        if (!blockExistsB[blockIndex]) continue;

        // get the block in B
        const MATRIX3 blockB = B.block(3 * k, 3 * j, 3,3);
        const MATRIX3 blockC = blockA * blockB;

        for (int y = 0; y < 3; y++)
          for (int x = 0; x < 3; x++)
          {
            REAL entry = blockC(x,y);
            TRIPLET triplet(3 * k + x, 3 * j + y, entry);

            triplets.push_back(triplet);
          }
      }
  }
  multiply.stop();
  
  TIMER assembly("Assembly");
  const int& totalRows = A._totalRows;
  SPARSE_MATRIX C(totalRows, totalRows);
  C.setFromTriplets(triplets.begin(), triplets.end());
  assembly.stop();
  return C;

  /*
  typedef Eigen::Triplet<REAL> TRIPLET;
  std::vector<TRIPLET> triplets;

  // it's square, right?
  assert(totalRows == A._totalCols);
  for (int i = 0; i < totalRows; i++)
    for (unsigned int x = 0; x < B.outerSize(); x++)
      for (SPARSE_MATRIX::InnerIterator it(B, x); it; ++it)
      {
        const int k = it.row();
        const int j = it.col();

        if (!A.exists(i,k)) continue;

        //C(i,j) += A.entry(i,k) * B(k,j);
        REAL entry = A.entry(i,k) * it.value();

        assert(i < totalRows);
        assert(j < totalRows);
        assert(i >= 0);
        assert(j >= 0);

        TRIPLET triplet(i,j,entry);
        triplets.push_back(triplet);
      }

  SPARSE_MATRIX C(totalRows, totalRows);
  C.setFromTriplets(triplets.begin(), triplets.end());

  return C;
  */
}
  /*
{
  TIMER functionTimer(__FUNCTION__);
  assert(A._totalRows == B.rows());
  assert(A._totalCols == B.cols());

  const int& totalRows = A._totalRows;

  typedef Eigen::Triplet<REAL> TRIPLET;
  std::vector<TRIPLET> triplets;

  // it's square, right?
  assert(totalRows == A._totalCols);
  for (int i = 0; i < totalRows; i++)
    for (unsigned int x = 0; x < B.outerSize(); x++)
      for (SPARSE_MATRIX::InnerIterator it(B, x); it; ++it)
      {
        const int k = it.row();
        const int j = it.col();

        if (!A.exists(i,k)) continue;

        //C(i,j) += A.entry(i,k) * B(k,j);
        REAL entry = A.entry(i,k) * it.value();

        assert(i < totalRows);
        assert(j < totalRows);
        assert(i >= 0);
        assert(j >= 0);

        TRIPLET triplet(i,j,entry);
        triplets.push_back(triplet);
      }

  SPARSE_MATRIX C(totalRows, totalRows);
  C.setFromTriplets(triplets.begin(), triplets.end());

  return C;
}
*/

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX3& BLOCK_DIAGONAL_MATRIX3::operator=(const BLOCK_DIAGONAL_MATRIX3& A)
{
  _totalRows = A._totalRows;
  _totalCols = A._totalCols;
  
  _totalBlockRows = A._totalBlockRows;
  _totalBlockCols = A._totalBlockCols;

  _blocks.resize(_totalBlockRows);
  for (int x = 0; x < _totalBlockRows; x++)
    _blocks[x] = A._blocks[x];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// set to identity matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_DIAGONAL_MATRIX3::setIdentity()
{
  for (int x = 0; x < _totalBlockRows; x++)
    _blocks[x].setIdentity();
}

