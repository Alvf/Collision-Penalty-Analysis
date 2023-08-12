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

#ifndef BLOCK_SPARSE_MATRIX3_H
#define BLOCK_SPARSE_MATRIX3_H

#include "SETTINGS.h"
#include "util/BLOCK_DIAGONAL_MATRIX3.h"
//#include <map>
#include <unordered_map>

class BLOCK_SPARSE_MATRIX3 
{
public:
  BLOCK_SPARSE_MATRIX3(const int& totalRows, const int& totalCols);
  BLOCK_SPARSE_MATRIX3(const SPARSE_MATRIX& A);

  // access a block matrix
  MATRIX3& operator()(const int row, const int col);

  const int rows() const { return _totalRows; };
  const int cols() const { return _totalCols; };

  // take the transpose
  BLOCK_SPARSE_MATRIX3 transpose() const;
 
  // scale on both the left and right
  BLOCK_SPARSE_MATRIX3 scaleLeftRight(const BLOCK_DIAGONAL_MATRIX3& A);

  // scale on both the left and right and then add B
  BLOCK_SPARSE_MATRIX3 scaleLeftRightAdd(const BLOCK_DIAGONAL_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B);

  friend BLOCK_SPARSE_MATRIX3 operator+(const BLOCK_SPARSE_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B);
  friend BLOCK_SPARSE_MATRIX3 operator*(const BLOCK_SPARSE_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B);
  friend BLOCK_SPARSE_MATRIX3 operator*(const BLOCK_DIAGONAL_MATRIX3& A, const BLOCK_SPARSE_MATRIX3& B);

  // UNTESTED
  //friend VECTOR operator*(const BLOCK_SPARSE_MATRIX3& A, const VECTOR& v);

  SPARSE_MATRIX sparseView() const;

  const bool blockExists(const int blockRow, const int blockCol) const;

  // wipe everything
  void setZero();

protected:
  int _totalRows;
  int _totalCols;
  int _totalBlockRows;
  int _totalBlockCols;

  int _totalBlocks;
  int _totalNonZeroBlocks;
  //std::map<int, MATRIX3> _blocks;
  //std::vector<std::map<int, MATRIX3> > _blocks;
  std::vector<std::unordered_map<int, MATRIX3> > _blocks;
  //std::map<int, bool> _blockExists;
  std::unordered_map<int, bool> _blockExists;
};

#endif
