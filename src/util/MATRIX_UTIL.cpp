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
#include "MATRIX_UTIL.h"
#include "TIMER.h"

#include <random>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

std::mt19937 gen(123);
std::uniform_real_distribution<REAL> dist(0.0, 1.0);

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// convert a MATRIX3 to a VECTOR9 in a consistent way
///////////////////////////////////////////////////////////////////////
VECTOR9 flatten(const MATRIX3& A)
{
  VECTOR9 column;

  unsigned int index = 0;
  for (unsigned int j = 0; j < A.cols(); j++)
    for (unsigned int i = 0; i < A.rows(); i++, index++)
      column[index] = A(i,j);

  return column;
}

///////////////////////////////////////////////////////////////////////
// convert a MATRIX3x2 to a VECTOR6 in a consistent way
///////////////////////////////////////////////////////////////////////
VECTOR6 flatten(const MATRIX3x2& A)
{
  VECTOR6 column;

  unsigned int index = 0;
  for (unsigned int j = 0; j < A.cols(); j++)
    for (unsigned int i = 0; i < A.rows(); i++, index++)
      column[index] = A(i,j);

  return column;
}

///////////////////////////////////////////////////////////////////////
// convert a VECTOR9 to a MATRIX3 in a consistent way
///////////////////////////////////////////////////////////////////////
MATRIX3 unflatten(const VECTOR9& v)
{
  MATRIX3 A;
  unsigned int index = 0;
  for (unsigned int j = 0; j < A.cols(); j++)
    for (unsigned int i = 0; i < A.rows(); i++, index++)
      A(i,j) = v[index];

  return A;
}

///////////////////////////////////////////////////////////////////////
// get the polar decomposition of matrix A = RS
///////////////////////////////////////////////////////////////////////
void polarDecomposition(const MATRIX3& A, MATRIX3& R, MATRIX3& S)
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(A, U, Sigma, V);

  R = U * V.transpose();
  S = V * Sigma.asDiagonal() * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// get the polar decomposition of matrix A = RS
///////////////////////////////////////////////////////////////////////
void polarDecomposition(const MATRIX3x2& A, MATRIX3x2& R, MATRIX2& S)
{
  MATRIX3x2 U;
  VECTOR2 sigma;
  MATRIX2 V;
  svd(A, U, sigma, V);

  R = U * V.transpose();
  S = V * sigma.asDiagonal() * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// rotation variant of the SVD where the reflections are loaded into
// Sigma and not U and V
///////////////////////////////////////////////////////////////////////
void svd_rv(const MATRIX3& F, MATRIX3& U, VECTOR3& Sigma, MATRIX3& V)
{
  const Eigen::JacobiSVD<MATRIX3,Eigen::NoQRPreconditioner> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  Sigma = svd.singularValues();

  MATRIX3 L = MATRIX3::Identity();
  L(2,2) = (U * V.transpose()).determinant();

  const REAL detU = U.determinant();
  const REAL detV = V.determinant();

  if (detU < 0.0 && detV > 0)
    U = U * L;
  if (detU > 0.0 && detV < 0.0)
    V = V * L;

  Sigma[2] = Sigma[2] * L(2,2);
}

///////////////////////////////////////////////////////////////////////
// svd of a 3x2, which has no rotation variant, because a reflection 
// just means a pi rotation
///////////////////////////////////////////////////////////////////////
void svd(const MATRIX3x2& F, MATRIX3x2& U, VECTOR2& Sigma, MATRIX2& V)
{
  using namespace Eigen;

  // casting here is a bit troubling
  //const Eigen::JacobiSVD<MATRIX3x2,Eigen::NoQRPreconditioner> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
  //const JacobiSVD<MATRIX3x2> svd(F, ComputeThinU | ComputeThinV);
  const JacobiSVD<MATRIX> svd(F, ComputeThinU | ComputeThinV);
  U = svd.matrixU();
  V = svd.matrixV();
  Sigma = svd.singularValues();
}

///////////////////////////////////////////////////////////////////////
// clamp small values directly to zero, mostly just for printing
///////////////////////////////////////////////////////////////////////
MATRIX clampSmalls(const MATRIX& A, const REAL delta)
{
  MATRIX result = A;
  for (unsigned int y = 0; y < A.cols(); y++)
    for (unsigned int x = 0; x < A.rows(); x++)
      result(x,y) = fabs(result(x,y)) < delta ? 0.0 : result(x,y);

  return result;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 9x9 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX9 clampEigenvalues(const MATRIX9& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX9> eigensolver(A);
  const MATRIX9 Q = eigensolver.eigenvectors();
  VECTOR9 values = eigensolver.eigenvalues();
  for (int x = 0; x < 9; x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX9 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 6x6 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX6 clampEigenvalues(const MATRIX6& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX6> eigensolver(A);
  const MATRIX6 Q = eigensolver.eigenvectors();
  VECTOR6 values = eigensolver.eigenvalues();
  for (int x = 0; x < 6; x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX6 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 3x3 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX3 clampEigenvalues(const MATRIX3& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX3> eigensolver(A);
  const MATRIX3 Q = eigensolver.eigenvectors();
  VECTOR3 values = eigensolver.eigenvalues();
  for (int x = 0; x < 3; x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX3 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 4x4 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX4 clampEigenvalues(const MATRIX4& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX4> eigensolver(A);
  const MATRIX4 Q = eigensolver.eigenvectors();
  VECTOR4 values = eigensolver.eigenvalues();
  for (int x = 0; x < 4; x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX4 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 11x11 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX11 clampEigenvalues(const MATRIX11& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX11> eigensolver(A);
  const MATRIX11 Q = eigensolver.eigenvectors();
  VECTOR11 values = eigensolver.eigenvalues();
  for (int x = 0; x < 11; x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX11 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 12x12 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX12 clampEigenvalues(const MATRIX12& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX12> eigensolver(A);
  const MATRIX12 Q = eigensolver.eigenvectors();
  VECTOR12 values = eigensolver.eigenvalues();
  for (int x = 0; x < 12; x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX12 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 12x12 to semi-negative-definite
///////////////////////////////////////////////////////////////////////
MATRIX12 clampEigenvaluesToSemiNegative(const MATRIX12& A)
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX12> eigensolver(A);
  const MATRIX12 Q = eigensolver.eigenvectors();
  VECTOR12 values = eigensolver.eigenvalues();
  for (int x = 0; x < 12; x++)
    values[x] = (values[x] < 0.0) ? values[x] : 0.0;
  MATRIX12 B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 2x2 matrix
///////////////////////////////////////////////////////////////////////
void eigensystem(const MATRIX2& A, MATRIX2& Q, VECTOR2& Lambda)
{
  Eigen::SelfAdjointEigenSolver<MATRIX2> eigensolver(A);
  Q = eigensolver.eigenvectors();
  Lambda = eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 3x3 matrix
///////////////////////////////////////////////////////////////////////
void eigensystem(const MATRIX3& A, MATRIX3& Q, VECTOR3& Lambda)
{
  Eigen::SelfAdjointEigenSolver<MATRIX3> eigensolver(A);
  Q = eigensolver.eigenvectors();
  Lambda = eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 6x6 matrix
///////////////////////////////////////////////////////////////////////
void eigensystem(const MATRIX6& A, MATRIX6& Q, VECTOR6& Lambda)
{
  Eigen::SelfAdjointEigenSolver<MATRIX6> eigensolver(A);
  Q = eigensolver.eigenvectors();
  Lambda = eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 9x9 matrix
///////////////////////////////////////////////////////////////////////
void eigensystem(const MATRIX9& A, MATRIX9& Q, VECTOR9& Lambda)
{
  Eigen::SelfAdjointEigenSolver<MATRIX9> eigensolver(A);
  Q = eigensolver.eigenvectors();
  Lambda = eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 12x12 matrix
///////////////////////////////////////////////////////////////////////
void eigensystem(const MATRIX12& A, MATRIX12& Q, VECTOR12& Lambda)
{
  Eigen::SelfAdjointEigenSolver<MATRIX12> eigensolver(A);
  Q = eigensolver.eigenvectors();
  Lambda = eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 3x3 matrix
///////////////////////////////////////////////////////////////////////
VECTOR3 eigenvalues(const MATRIX3& A)
{
  Eigen::SelfAdjointEigenSolver<MATRIX3> eigensolver(A);
  return eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 6x6 matrix
///////////////////////////////////////////////////////////////////////
VECTOR6 eigenvalues(const MATRIX6& A)
{
  Eigen::SelfAdjointEigenSolver<MATRIX6> eigensolver(A);
  return eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 9x9 matrix
///////////////////////////////////////////////////////////////////////
VECTOR9 eigenvalues(const MATRIX9& A)
{
  Eigen::SelfAdjointEigenSolver<MATRIX9> eigensolver(A);
  return eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of a 12x12 matrix
///////////////////////////////////////////////////////////////////////
VECTOR12 eigenvalues(const MATRIX12& A)
{
  Eigen::SelfAdjointEigenSolver<MATRIX12> eigensolver(A);
  return eigensolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////
// get the eigensystem of an arbitrary matrix
///////////////////////////////////////////////////////////////////////
VECTOR eigenvalues(const MATRIX& A)
{
  Eigen::SelfAdjointEigenSolver<MATRIX> eigensolver(A);
  return eigensolver.eigenvalues();
} 

///////////////////////////////////////////////////////////////////////
// get the eigensystem of an arbitrary sparse matrix
///////////////////////////////////////////////////////////////////////
VECTOR eigenvalues(const SPARSE_MATRIX& A)
{ 
  return eigenvalues(MATRIX(A)); 
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL largestEigenvalue(const SPARSE_MATRIX& A)
{
  using namespace Spectra;
  
  // Construct matrix operation object using the wrapper class SparseGenMatProd
  SparseSymMatProd<REAL> op(A);

  // Construct eigen solver object, requesting the largest three eigenvalues
  SymEigsSolver<REAL, LARGEST_MAGN, SparseSymMatProd<REAL> > eigs(&op, 3, 6);

  // Initialize and compute
  eigs.init();
  eigs.compute();

  // Retrieve results
  VECTOR evalues;
  evalues = eigs.eigenvalues();

  return evalues[0];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL smallestEigenvalue(const SPARSE_MATRIX& A)
{
  using namespace Spectra;

  // Construct matrix operation object using the wrapper class SparseGenMatProd
  SparseSymShiftSolve<REAL> op(A);

  // Construct eigen solver object, requesting the largest three eigenvalues
  SymEigsShiftSolver<REAL, SMALLEST_MAGN, SparseSymShiftSolve<REAL> > eigs(&op, 3, 6, 0.0);

  // Initialize and compute
  eigs.init();
  eigs.compute();

  // Retrieve results
  VECTOR evalues;
  evalues = eigs.eigenvalues();

  return evalues[0];
}

///////////////////////////////////////////////////////////////////////
// Let's make some random deformation gradients
///////////////////////////////////////////////////////////////////////
MATRIX3 randomMatrix3(const REAL scaling)
{
  MATRIX3 F;
  F.setIdentity();

  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
    {
      F(x,y) = dist(gen);

      // randomize the sign
      if (dist(gen) < 0.5)
        F(x,y) *= -1.0;
    }
  F *= scaling;
  return F;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random deformation gradients
///////////////////////////////////////////////////////////////////////
MATRIX3x2 randomMatrix3x2(const REAL scaling)
{
  MATRIX3x2 F;
  F.setIdentity();

  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 3; x++)
    {
      F(x,y) = dist(gen);

      // randomize the sign
      if (dist(gen) < 0.5)
        F(x,y) *= -1.0;
    }
  F *= scaling;
  return F;
}

MATRIX4 randomMatrix4(const REAL scaling)
{
  MATRIX4 F;
  F.setIdentity();

  for (int y = 0; y < 4; y++)
    for (int x = 0; x < 4; x++)
    {
      F(x,y) = dist(gen);

      // randomize the sign
      if (dist(gen) < 0.5)
        F(x,y) *= -1.0;
    }
  F *= scaling;
  return F;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random positive-definite deformation gradients
///////////////////////////////////////////////////////////////////////
MATRIX3 randomPositiveDefiniteMatrix3(const REAL scaling)
{
  MATRIX3 U = randomRotation();
  MATRIX3 V = randomRotation();

  MATRIX3 Sigma;
  Sigma.setIdentity();

  for (int x = 0; x < 3; x++)
    Sigma(x,x) = dist(gen);
  Sigma *= scaling;
  return U * Sigma * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// Let's make a random barycentric coordinate
///////////////////////////////////////////////////////////////////////
VECTOR2 randomBarycentric()
{
  VECTOR2 v;
  v[0] = dist(gen);
  v[1] = 1.0 - v[0];
  return v;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random directions
///////////////////////////////////////////////////////////////////////
VECTOR3 randomVector3(const REAL scaling)
{
  VECTOR3 v;
  for (int x = 0; x < 3; x++)
  {
    v[x] = dist(gen);

    // randomize the sign
    if (dist(gen) < 0.5)
      v[x] *= -1.0;
  }
  v *= scaling;
  return v;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random directions
///////////////////////////////////////////////////////////////////////
VECTOR12 randomVector12(const REAL scaling)
{
  VECTOR12 v;
  for (int x = 0; x < 12; x++)
  {
    v[x] = dist(gen);

    // randomize the sign
    if (dist(gen) < 0.5)
      v[x] *= -1.0;
  }
  v *= scaling;
  return v;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random rotations
///////////////////////////////////////////////////////////////////////
MATRIX3 randomRotation()
{
  const REAL angle = dist(gen) * 2.0 * M_PI;
  const VECTOR3 axis = randomVector3().normalized();
  MATRIX3 R;

  // assumes that REAL is a double
  R = Eigen::AngleAxisd(angle, axis);
  return R;
}

///////////////////////////////////////////////////////////////////////
// Matrix double-contraction
///////////////////////////////////////////////////////////////////////
REAL ddot(const MATRIX3& A, const MATRIX3& B)
{
  REAL result = 0;
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      result += A(x,y) * B(x,y);

  return result;
}

///////////////////////////////////////////////////////////////////////
// rotation gradient, w.r.t. deformation gradient F
// \frac{\partial R}{\partial F}
///////////////////////////////////////////////////////////////////////
MATRIX9 rotationGradient(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V)
{
  const REAL& sx = Sigma[0];
  const REAL& sy = Sigma[1];
  const REAL& sz = Sigma[2];

  // create the pseudo-twist vectors
  MATRIX3 twistx, twisty, twistz; 
  twistx <<  0, 0, 0,
             0, 0, -1,
             0, 1, 0;
  twisty <<  0, 0, 1,
             0, 0, 0,
            -1, 0, 0;
  twistz <<  0, 1, 0,
            -1, 0, 0,
             0, 0, 0;

  // compute the eigenvectors
  const REAL front = 1.0 / sqrt(2.0);
  const MATRIX3 Qx = front * U * twistx * V.transpose();
  const MATRIX3 Qy = front * U * twisty * V.transpose();
  const MATRIX3 Qz = front * U * twistz * V.transpose();

  // flatten them out to vectors
  const VECTOR9 qx = flatten(Qx);
  const VECTOR9 qy = flatten(Qy);
  const VECTOR9 qz = flatten(Qz);

  // compute the eigenvectors
  const REAL lambdax = 2.0 / (sy + sz);
  const REAL lambday = 2.0 / (sx + sz);
  const REAL lambdaz = 2.0 / (sx + sy);

  MATRIX9 gradient;

  gradient  = lambdax * (qx * qx.transpose());
  gradient += lambday * (qy * qy.transpose());
  gradient += lambdaz * (qz * qz.transpose());
  return gradient;
}

///////////////////////////////////////////////////////////////////////
// time derivative of rotation
///////////////////////////////////////////////////////////////////////
MATRIX3 rotationDot(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V, const MATRIX3& Fdot)
{
  MATRIX9 DRDF = rotationGradient(U, Sigma, V);
  VECTOR9 fdot = flatten(Fdot);
  VECTOR9 rdot = DRDF * fdot;

  return unflatten(rdot);
}

///////////////////////////////////////////////////////////////////////
// eigenvectors 0-2 are the twist modes
// eigenvectors 3-5 are the flip modes
///////////////////////////////////////////////////////////////////////
void buildTwistAndFlipEigenvectors(const MATRIX3& U, const MATRIX3& V, MATRIX9& Q)
{
  // create the twist matrices
  MATRIX3 T0, T1, T2; 
  T0 <<  0, 0, 0,
         0, 0, -1,
         0, 1, 0;   // x-twist
  T1 <<  0, 0, 1,
         0, 0, 0,
        -1, 0, 0;   // y-twist
  T2 <<  0, 1, 0,
        -1, 0, 0,
         0, 0, 0;   // z-twist

  const MATRIX3 Q0 = (1.0 / sqrt(2.0)) * (U * T0 * V.transpose());
  const MATRIX3 Q1 = (1.0 / sqrt(2.0)) * (U * T1 * V.transpose());
  const MATRIX3 Q2 = (1.0 / sqrt(2.0)) * (U * T2 * V.transpose());

  // create the flip matrices
  MATRIX3 L0, L1, L2; 
  L0 <<  0, 0, 0,
         0, 0, 1,
         0, 1, 0;   // x-flip
  L1 <<  0, 0, 1,
         0, 0, 0,
         1, 0, 0;   // y-flip
  L2 <<  0, 1, 0,
         1, 0, 0,
         0, 0, 0;   // z-flip

  const MATRIX3 Q3 = (1.0 / sqrt(2.0)) * (U * L0 * V.transpose());
  const MATRIX3 Q4 = (1.0 / sqrt(2.0)) * (U * L1 * V.transpose());
  const MATRIX3 Q5 = (1.0 / sqrt(2.0)) * (U * L2 * V.transpose());

  Q.col(0) = flatten(Q0);
  Q.col(1) = flatten(Q1);
  Q.col(2) = flatten(Q2);
  Q.col(3) = flatten(Q3);
  Q.col(4) = flatten(Q4);
  Q.col(5) = flatten(Q5);
}

///////////////////////////////////////////////////////////////////////
// eigenvectors 6-8 are the scaling modes, non-jackpot version
///////////////////////////////////////////////////////////////////////
void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& Q,
                              const MATRIX3& V, MATRIX9& Q9)
{
  const VECTOR3 q0 = Q.col(0);
  const VECTOR3 q1 = Q.col(1);
  const VECTOR3 q2 = Q.col(2);
  
  const MATRIX3 Q0 = U * q0.asDiagonal() * V.transpose();
  const MATRIX3 Q1 = U * q1.asDiagonal() * V.transpose();
  const MATRIX3 Q2 = U * q2.asDiagonal() * V.transpose();

  Q9.col(6) = flatten(Q0);
  Q9.col(7) = flatten(Q1);
  Q9.col(8) = flatten(Q2);
}

///////////////////////////////////////////////////////////////////////
// eigenvectors 6-8 are the scaling modes, jackpot version
///////////////////////////////////////////////////////////////////////
void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& V, MATRIX9& Q9)
{
  VECTOR3 x(1,0,0);
  VECTOR3 y(0,1,0);
  VECTOR3 z(0,0,1);
  
  const MATRIX3 Q0 = U * x.asDiagonal() * V.transpose();
  const MATRIX3 Q1 = U * y.asDiagonal() * V.transpose();
  const MATRIX3 Q2 = U * z.asDiagonal() * V.transpose();

  Q9.col(6) = flatten(Q0);
  Q9.col(7) = flatten(Q1);
  Q9.col(8) = flatten(Q2);
}

///////////////////////////////////////////////////////////////////////
// get the Kronecker product of matrix with respect to 3x3 identity,
// used a lot in anisotropic materials
//
// in Matlab, 
//
// I = eye(3,3);
// H = [A(1,1) * I A(1,2) * I A(1,3) * I;
//      A(2,1) * I A(2,2) * I A(2,3) * I;
//      A(3,1) * I A(3,2) * I A(3,3) * I];
//
// or more succinctly: kron(A,eye(3,3)) 
///////////////////////////////////////////////////////////////////////
MATRIX9 kronIdentity(const MATRIX3& A)
{
  MATRIX9 H = MATRIX9::Zero();
  H.block<3,3>(0,0) = MATRIX3::Identity() * A(0,0);
  H.block<3,3>(3,3) = MATRIX3::Identity() * A(1,1);
  H.block<3,3>(6,6) = MATRIX3::Identity() * A(2,2);

  H.block<3,3>(3,0) = H.block<3,3>(0,3) = MATRIX3::Identity() * A(1,0);
  H.block<3,3>(6,0) = H.block<3,3>(0,6) = MATRIX3::Identity() * A(2,0);
  H.block<3,3>(6,3) = H.block<3,3>(3,6) = MATRIX3::Identity() * A(1,2);

  return H;
}

///////////////////////////////////////////////////////////////////////
// get the Kronecker product of matrix with respect to 3x3 identity,
// used a lot in anisotropic materials
//
// in Matlab, 
//
// I = eye(3,3);
// H = [A(1,1) * I A(1,2) * I 
//      A(2,1) * I A(2,2) * I]; 
//
// or more succinctly: kron(A,eye(3,3)) 
///////////////////////////////////////////////////////////////////////
MATRIX6 kronIdentity(const MATRIX2& A)
{
  MATRIX6 H = MATRIX6::Zero();
  H.block<3,3>(0,0) = MATRIX3::Identity() * A(0,0);
  H.block<3,3>(3,3) = MATRIX3::Identity() * A(1,1);
  H.block<3,3>(3,0) = H.block<3,3>(0,3) = MATRIX3::Identity() * A(1,0);
  return H;
}

///////////////////////////////////////////////////////////////////////
// Tensor invariants, 3D volumes
///////////////////////////////////////////////////////////////////////
REAL invariant2(const MATRIX3& F)
{
  return ddot(F,F);
}
REAL invariant2(const VECTOR3& Sigma)
{
  return Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1] + Sigma[2] * Sigma[2];
}
REAL invariant3(const MATRIX3& F)
{
  return F.determinant();
}
REAL invariant3(const VECTOR3& Sigma)
{
  return Sigma[0] * Sigma[1] * Sigma[2];
}
REAL invariant4(const MATRIX3& F, const VECTOR3& a)
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  const MATRIX3 S = V * Sigma.asDiagonal() * V.transpose();
  return (S * a).dot(a);
}
REAL invariant5(const MATRIX3& F, const VECTOR3& a)
{
  return (F * a).squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// Tensor invariants, 2D membranes
///////////////////////////////////////////////////////////////////////
REAL invariant1(const MATRIX3x2& F)
{
  MATRIX3x2 R;
  MATRIX2 S;
  polarDecomposition(F, R, S);
  return S.trace();
}
REAL invariant2(const MATRIX3x2& F)
{
  return F.squaredNorm();
}
REAL invariant3(const MATRIX3x2& F)
{
  MATRIX3x2 U;
  MATRIX2 V;
  VECTOR2 Sigma;
  svd(F, U, Sigma, V);
  return Sigma[0] * Sigma[1];
}

//////////////////////////////////////////////////////////////////////////////
// Eqn. 19 from Section 4.2 in "Stable Neo-Hookean Flesh Simulation"
//////////////////////////////////////////////////////////////////////////////
MATRIX3 partialJpartialF(const MATRIX3& F)
{
  MATRIX3 pJpF;
  pJpF.col(0) = F.col(1).cross(F.col(2));
  pJpF.col(1) = F.col(2).cross(F.col(0));
  pJpF.col(2) = F.col(0).cross(F.col(1));
  return pJpF;
}

//////////////////////////////////////////////////////////////////////////////
// Eqn. 29 from Section 4.5 in "Stable Neo-Hookean Flesh Simulation",
// with a scaling factor added
//////////////////////////////////////////////////////////////////////////////
MATRIX3 crossProduct(const MATRIX3& F, const int i)
{
  return (MATRIX(3,3) <<       0, -F(2,i),  F(1,i),
                          F(2,i),       0, -F(0,i),
                         -F(1,i),  F(0,i),       0).finished();
}

//////////////////////////////////////////////////////////////////////////////
// cross product matrix
//////////////////////////////////////////////////////////////////////////////
MATRIX3 crossProduct(const VECTOR3& x)
{
  return (MATRIX(3,3) <<     0, -x[2],  x[1],
                          x[2],     0, -x[0],
                         -x[1],  x[0],     0).finished();
}

//////////////////////////////////////////////////////////////////////////////
// 3rd order tensor derivative of deformation gradient F with respect 
// to itself
//////////////////////////////////////////////////////////////////////////////
void partialFpartialF(const int i, const int j, MATRIX3& pFpF)
{
  pFpF.setZero();
  pFpF(i,j) = 1;
}
  
//////////////////////////////////////////////////////////////////////////////
// 3D axis-angle rotation matrix
// angle is in radians
//////////////////////////////////////////////////////////////////////////////
MATRIX3 rotationMatrix(const VECTOR3& axis, const REAL& angle)
{
  if (axis.norm() < 1e-8) return MATRIX3::Identity();

  const VECTOR3& n = axis.normalized();

  return Eigen::AngleAxisd(angle, axis).toRotationMatrix();
}

//////////////////////////////////////////////////////////////////////////////
// reset the seed in the randomization variable
//////////////////////////////////////////////////////////////////////////////
void reseed(const int seedValue)
{
  gen.seed(seedValue);
}

} // HOBAK
