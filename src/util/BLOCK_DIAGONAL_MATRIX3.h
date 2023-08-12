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

#ifndef BLOCK_DIAGONAL_MATRIX3_H
#define BLOCK_DIAGONAL_MATRIX3_H

#include "SETTINGS.h"
#include <vector>

class BLOCK_DIAGONAL_MATRIX3 
{
public:
  BLOCK_DIAGONAL_MATRIX3();
  BLOCK_DIAGONAL_MATRIX3(const int& totalRows, const int& totalCols);

  // access a block matrix
  MATRIX3& operator()(const int blockRow, const int blockCol);
  const MATRIX3& operator()(const int blockRow, const int blockCol) const;
  REAL& entry(const int row, const int col);
  const REAL entry(const int row, const int col) const;

  BLOCK_DIAGONAL_MATRIX3& operator=(const BLOCK_DIAGONAL_MATRIX3& A);

  // is there even an entry to retrieve?
  const bool exists(const int row, const int col) const;

  // set to identity matrix
  void setIdentity();

  /*
  const int& totalBlockRows() const { return _totalBlockRows; };
  const int& totalBlockCols() const { return _totalBlockCols; };
  */
  const int& rows() const { return _totalRows; };
  const int& cols() const { return _totalCols; };

  // return a transposed version
  BLOCK_DIAGONAL_MATRIX3 transpose() const;

  //friend SPARSE_MATRIX operator*(const SPARSE_MATRIX& A, const BLOCK_DIAGONAL_MATRIX3& B);
  friend SPARSE_MATRIX operator*(const BLOCK_DIAGONAL_MATRIX3& A, const SPARSE_MATRIX& B);
  friend SPARSE_MATRIX operator+(const SPARSE_MATRIX& A, const BLOCK_DIAGONAL_MATRIX3& B);
  friend BLOCK_DIAGONAL_MATRIX3 operator-(const BLOCK_DIAGONAL_MATRIX3& A, const BLOCK_DIAGONAL_MATRIX3& B);
  friend VECTOR operator*(const BLOCK_DIAGONAL_MATRIX3& A, const VECTOR& v);

protected:
  int _totalRows;
  int _totalCols;
  int _totalBlockRows;
  int _totalBlockCols;

  std::vector<MATRIX3> _blocks;
};

#endif
