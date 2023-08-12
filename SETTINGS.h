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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef double REAL;

typedef Eigen::Matrix<REAL, 2,  2>  MATRIX2;
typedef Eigen::Matrix<REAL, 3,  2>  MATRIX3x2;
typedef Eigen::Matrix<REAL, 3,  3>  MATRIX3;
typedef Eigen::Matrix<REAL, 3,  6>  MATRIX3x6;
typedef Eigen::Matrix<REAL, 3,  12> MATRIX3x12;
typedef Eigen::Matrix<REAL, 4,  4>  MATRIX4;
typedef Eigen::Matrix<REAL, 6,  6>  MATRIX6;
typedef Eigen::Matrix<REAL, 6,  4>  MATRIX6x4;
typedef Eigen::Matrix<REAL, 6,  9>  MATRIX6x9;
typedef Eigen::Matrix<REAL, 6,  11>  MATRIX6x11;
typedef Eigen::Matrix<REAL, 6,  12> MATRIX6x12;
typedef Eigen::Matrix<REAL, 9,  9>  MATRIX9;
typedef Eigen::Matrix<REAL, 9,  4>  MATRIX9x4;
typedef Eigen::Matrix<REAL, 9,  12> MATRIX9x12;
typedef Eigen::Matrix<REAL, 11,  2>  MATRIX11x2;
typedef Eigen::Matrix<REAL, 11, 11>  MATRIX11;
typedef Eigen::Matrix<REAL, 12, 12> MATRIX12;
typedef Eigen::Matrix<REAL, 2,  1>  VECTOR2;
typedef Eigen::Matrix<REAL, 3,  1>  VECTOR3;
typedef Eigen::Matrix<REAL, 4,  1>  VECTOR4;
typedef Eigen::Matrix<REAL, 6,  1>  VECTOR6;
typedef Eigen::Matrix<REAL, 9,  1>  VECTOR9;
typedef Eigen::Matrix<REAL, 11,  1>  VECTOR11;
typedef Eigen::Matrix<REAL, 12, 1>  VECTOR12;

typedef Eigen::Matrix<int, 2, 1> VECTOR2I;
typedef Eigen::Matrix<int, 3, 1> VECTOR3I;
typedef Eigen::Matrix<int, 4, 1> VECTOR4I;
typedef Eigen::Matrix<int, 5, 1> VECTOR5I;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VECTORI;

typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> VECTOR;
typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> MATRIX;
//typedef Eigen::SparseMatrix<REAL, Eigen::ColMajor> SPARSE_MATRIX;
typedef Eigen::SparseMatrix<REAL, Eigen::RowMajor> SPARSE_MATRIX;

// print out more debug information than usual?
#define VERY_VERBOSE 0

// #define CHEAT_POKE_THROUGH 1

// enable debugging traps?
#define ENABLE_DEBUG_TRAPS 0

#endif
