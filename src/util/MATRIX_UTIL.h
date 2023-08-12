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
#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H

#include "SETTINGS.h"

namespace HOBAK {

  // convert a MATRIX3 to a VECTOR9 in a consistent way
  VECTOR9 flatten(const MATRIX3& A);

  // convert a MATRIX3x2 to a VECTOR6 in a consistent way
  VECTOR6 flatten(const MATRIX3x2& A);
  
  // convert a VECTOR9 to a MATRIX3 in a consistent way
  MATRIX3 unflatten(const VECTOR9& v);

  // rotation variant of the SVD where the reflections are loaded into
  // Sigma and not U and V
  void svd_rv(const MATRIX3& F, MATRIX3& U, VECTOR3& Sigma, MATRIX3& V);

  // svd of a 3x2, which has no rotation variant, because a reflection just
  // means a pi rotation
  void svd(const MATRIX3x2& F, MATRIX3x2& U, VECTOR2& Sigma, MATRIX2& V);

  // get the polar decomposition of matrix A = RS
  void polarDecomposition(const MATRIX3& A, MATRIX3& R, MATRIX3& S);
  void polarDecomposition(const MATRIX3x2& A, MATRIX3x2& R, MATRIX2& S);

  // clamp the eigenvalues to semi-positive-definite
  MATRIX9 clampEigenvalues(const MATRIX9& A);
  MATRIX6 clampEigenvalues(const MATRIX6& A);
  MATRIX3 clampEigenvalues(const MATRIX3& A);
  MATRIX4 clampEigenvalues(const MATRIX4& A);
  MATRIX11 clampEigenvalues(const MATRIX11& A);
  MATRIX12 clampEigenvalues(const MATRIX12& A);
  MATRIX12 clampEigenvaluesToSemiNegative(const MATRIX12& A);

  // get the eigensystem of various matrices
  void eigensystem(const MATRIX2& A, MATRIX2& Q, VECTOR2& Lambda);
  void eigensystem(const MATRIX3& A, MATRIX3& Q, VECTOR3& Lambda);
  void eigensystem(const MATRIX6& A, MATRIX6& Q, VECTOR6& Lambda);
  void eigensystem(const MATRIX9& A, MATRIX9& Q, VECTOR9& Lambda);
  void eigensystem(const MATRIX12& A, MATRIX12& Q, VECTOR12& Lambda);
  VECTOR3 eigenvalues(const MATRIX3& A);
  VECTOR6 eigenvalues(const MATRIX6& A);
  VECTOR9 eigenvalues(const MATRIX9& A);
  VECTOR12 eigenvalues(const MATRIX12& A);
  VECTOR eigenvalues(const MATRIX& A);
  VECTOR eigenvalues(const SPARSE_MATRIX& A);

  // use Spectra to get the largest eig of a big matrix
  REAL largestEigenvalue(const SPARSE_MATRIX& A);
  
  // use Spectra to get the smallest eig of a big matrix
  REAL smallestEigenvalue(const SPARSE_MATRIX& A);

  // Let's make some random deformation gradients
  MATRIX3 randomMatrix3(const REAL scaling = 3.0);
  MATRIX3x2 randomMatrix3x2(const REAL scaling = 3.0);
  MATRIX4 randomMatrix4(const REAL scaling = 3.0);

  // Let's make some random positive-definite deformation gradients
  MATRIX3 randomPositiveDefiniteMatrix3(const REAL scaling = 3.0);

  // Let's make some random directions
  VECTOR3 randomVector3(const REAL scaling = 3.0);

  // Let's make some random directions
  VECTOR12 randomVector12(const REAL scaling = 10.0);

  // Let's make some random rotations
  MATRIX3 randomRotation();

  // Let's make a random barycentric coordinate
  VECTOR2 randomBarycentric();

  // Matrix double-contraction
  REAL ddot(const MATRIX3& A, const MATRIX3& B);

  // eigenvectors 0-2 are the twist modes
  // eigenvectors 3-5 are the flip modes
  void buildTwistAndFlipEigenvectors(const MATRIX3& U, const MATRIX3& V, 
                                     MATRIX9& Q);

  // eigenvectors 6-8 are the scaling modes, jackpot version
  void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& V, 
                                MATRIX9& Q9);

  // eigenvectors 6-8 are the scaling modes, non-jackpot version
  void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& Q,
                                const MATRIX3& V, MATRIX9& Q9);

  // get the Kronecker product of matrix with respect to 3x3 identity,
  // used a lot in anisotropic materials
  MATRIX9 kronIdentity(const MATRIX3& A);

  // get the Kronecker product of matrix with respect to 3x3 identity,
  // used with strand stretching
  MATRIX6 kronIdentity(const MATRIX2& A);

  // Tensor invariants, 3D volumes
  REAL invariant1(const MATRIX3& F);
  REAL invariant2(const MATRIX3& F);
  REAL invariant3(const MATRIX3& F);
  REAL invariant4(const MATRIX3& F, const VECTOR3& a);
  REAL invariant5(const MATRIX3& F, const VECTOR3& a);

  REAL invariant2(const VECTOR3& Sigma);
  REAL invariant3(const VECTOR3& Sigma);

  // Tensor invariants, 2D membranes
  REAL invariant1(const MATRIX3x2& F);
  REAL invariant2(const MATRIX3x2& F);
  REAL invariant3(const MATRIX3x2& F);

  // 3D axis-angle rotation matrix
  // angle is in radians
  MATRIX3 rotationMatrix(const VECTOR3& axis, const REAL& angle);

  // rotation gradient, w.r.t. deformation gradient F
  // \frac{\partial R}{\partial F}
  MATRIX9 rotationGradient(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V);

  // time derivative of rotation
  // \frac{\partial R}{\partial t}
  MATRIX3 rotationDot(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V, const MATRIX3& Fdot);

  // Eqn. 19 from Section 4.2 in "Stable Neo-Hookean Flesh Simulation"
  MATRIX3 partialJpartialF(const MATRIX3& F);

  // Eqn. 29 from Section 4.5 in "Stable Neo-Hookean Flesh Simulation"
  MATRIX3 crossProduct(const MATRIX3& F, const int col);
  
  // cross product matrix
  MATRIX3 crossProduct(const VECTOR3& x);

  // 3rd order tensor derivative of deformation gradient F with respect to itself
  void partialFpartialF(const int i, const int j, MATRIX3& pFpF);

  // clamp small values directly to zero, mostly just for printing
  MATRIX clampSmalls(const MATRIX& A, const REAL delta = 1e-7);
  
  // reset the seed in the randomization variable
  void reseed(const int seedValue);
}

#endif
