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
///////////////////////////////////////////////////////////////////////////////////////////////////
// This is a file in the HOBAK library
//
// July 11, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "QUADRATIC_UNIT_BENDING.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace STRAND {

QUADRATIC_UNIT_BENDING::QUADRATIC_UNIT_BENDING(const REAL& mu) :
  QUADRATIC_BENDING(mu)
{
}

QUADRATIC_UNIT_BENDING::QUADRATIC_UNIT_BENDING(const REAL& mu, const REAL& theta0) :
  QUADRATIC_BENDING(mu,theta0)
{
  //cout << " Bending theta0: " << theta0 << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string QUADRATIC_UNIT_BENDING::name() const
{ 
  return std::string("QUADRATIC_UNIT_BENDING"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL QUADRATIC_UNIT_BENDING::psi(const MATRIX3x2& E, const REAL& theta0) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const REAL diff = theta - theta0;

  return 0.5 * _mu * diff * diff;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL QUADRATIC_UNIT_BENDING::psi(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const REAL diff = (inverted) ? 2.0 * M_PI - theta - theta0 : theta - theta0;

  return 0.5 * _mu * diff * diff;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
static REAL clamp(const REAL& input, const REAL& bottom, const REAL& top)
{
  REAL result = (input > bottom) ? input : bottom;
  return result < top ? result : top;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 QUADRATIC_UNIT_BENDING::PK1(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  const REAL theta = acos(cosTheta);

  return _mu * (theta - theta0) * gradient(E);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 QUADRATIC_UNIT_BENDING::PK1(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  const REAL theta = acos(cosTheta);

  MATRIX3x2 Pfinal;
  if (!inverted)
    Pfinal = _mu * (theta - theta0) * gradient(E);
  else
    Pfinal = -_mu * (2.0 * M_PI - theta - theta0) * gradient(E);

  return Pfinal;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_UNIT_BENDING::hessian(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  const REAL cosTheta = t0.dot(t1);
  const REAL theta = acos(cosTheta);
  const REAL crossNorm = (t0.cross(t1)).norm();

  if (fabs(crossNorm) < 1e-8 || fabs(theta - theta0) < 1e-8)
  {
    // it will just be the outer product term, which is guaranteed semi-positive-definite
    VECTOR6 g = flatten(gradient(E));
    return _mu * g * g.transpose();
  }

  const VECTOR3 crossed = t0.cross(t1).normalized();
  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);

  VECTOR6 q2theta, q3theta, q4theta, q5theta;
  q2theta.setZero(); q3theta.setZero();
  q4theta.setZero(); q5theta.setZero();

  q2theta.block<3,1>(0,0) = t0 + e0cross;
  q3theta.block<3,1>(0,0) = t0 - e0cross;
  q4theta.block<3,1>(3,0) = t1 - e1cross;
  q5theta.block<3,1>(3,0) = t1 + e1cross;
  q2theta.normalize(); q3theta.normalize();
  q4theta.normalize(); q5theta.normalize();

  // segments calls bomb for some reason
  VECTOR6 q0, q1;
  q0.block<3,1>(0,0) = crossed;
  q0.block<3,1>(3,0) = crossed;
  q1.block<3,1>(0,0) = crossed;
  q1.block<3,1>(3,0) = -crossed;

  MATRIX3x2 T;
  T.col(0) = t0;
  T.col(1) = t1;
  const VECTOR6 p = flatten(gradient(T));

  const REAL dTheta = _mu * (theta - theta0);
  const REAL ddTheta = _mu;

  // build scaled version from eigendecomposition
  const REAL alpha = ddTheta * p.dot(p);

  VECTOR6 lambda;
  lambda[0] = dTheta * (cosTheta - 1.0) / crossNorm;
  lambda[1] = dTheta * (cosTheta + 1.0) / crossNorm;
  const REAL radical = sqrt(alpha * alpha + 4 * dTheta * dTheta);
  lambda[2] = (alpha + radical) * 0.5;
  lambda[3] = (alpha - radical) * 0.5;
  lambda[4] = dTheta;
  lambda[5] = -dTheta;

  VECTOR6 q2, q3, q4, q5;
  const REAL denom0 = (fabs(dTheta + lambda[2]) > 0.0) ? dTheta + lambda[2] : 1.0;
  const REAL denom1 = (fabs(dTheta - lambda[2]) > 0.0) ? dTheta - lambda[2] : 1.0;
  const REAL denom2 = (fabs(dTheta + lambda[3]) > 0.0) ? dTheta + lambda[3] : 1.0;
  const REAL denom3 = (fabs(dTheta - lambda[3]) > 0.0) ? dTheta - lambda[3] : 1.0;
  q2 = (1.0 / denom0) * (q2theta + q4theta) + (1.0 / denom1) * (q3theta + q5theta);
  q3 = (1.0 / denom2) * (q2theta + q4theta) + (1.0 / denom3) * (q3theta + q5theta);
  q4 = 0.5 * (q3theta - q5theta);
  q5 = 0.5 * (q2theta - q4theta);

  q0.normalize(); q1.normalize(); q2.normalize();
  q3.normalize(); q4.normalize(); q5.normalize();

  MATRIX6 S;
  S.setZero();
  S(0,0) = S(1,1) = S(2,2) = 1.0 / e0.norm();
  S(3,3) = S(4,4) = S(5,5) = 1.0 / e1.norm();

  //S = [eye(3,3) / e0norm zeros(3,3); zeros(3,3) eye(3,3) / e1norm];
  q0 = S * q0; q1 = S * q1; q2 = S * q2;
  q3 = S * q3; q4 = S * q4; q5 = S * q5;

  q0.normalize(); q1.normalize(); q2.normalize();
  q3.normalize(); q4.normalize(); q5.normalize();

  MATRIX6 Q;
  Q.col(0) = q0; Q.col(1) = q1; Q.col(2) = q2;
  Q.col(3) = q3; Q.col(4) = q4; Q.col(5) = q5;
  
  const REAL scale = 0.5 * (1 / e0.dot(e0) + 1 / e1.dot(e1)); 
  lambda = lambda * scale;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();
  return H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_UNIT_BENDING::hessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  const REAL cosTheta = t0.dot(t1);
  const REAL crossNorm = (t0.cross(t1)).norm();
  REAL theta = acos(cosTheta);

  if (inverted)
    theta = 2 * M_PI - theta;

  if (fabs(crossNorm) < 1e-8 || fabs(theta - theta0) < 1e-8)
  {
    // it will just be the outer product term, which is guaranteed semi-positive-definite
    VECTOR6 g = flatten(gradient(E));
    return _mu * g * g.transpose();
  }

  const VECTOR3 crossed = t0.cross(t1).normalized();
  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);

  VECTOR6 q2theta, q3theta, q4theta, q5theta;
  q2theta.setZero(); q3theta.setZero();
  q4theta.setZero(); q5theta.setZero();

  q2theta.block<3,1>(0,0) = t0 + e0cross;
  q3theta.block<3,1>(0,0) = t0 - e0cross;
  q4theta.block<3,1>(3,0) = t1 - e1cross;
  q5theta.block<3,1>(3,0) = t1 + e1cross;
  q2theta.normalize(); q3theta.normalize();
  q4theta.normalize(); q5theta.normalize();

  // segments calls bomb for some reason
  VECTOR6 q0, q1;
  q0.block<3,1>(0,0) = crossed;
  q0.block<3,1>(3,0) = crossed;
  q1.block<3,1>(0,0) = crossed;
  q1.block<3,1>(3,0) = -crossed;

  MATRIX3x2 T;
  T.col(0) = t0;
  T.col(1) = t1;
  const VECTOR6 p = flatten(gradient(T));

  REAL dTheta = _mu * (theta - theta0);
  REAL ddTheta = _mu;
  if (inverted)
    dTheta *= -1;

  // build scaled version from eigendecomposition
  const REAL alpha = ddTheta * p.dot(p);

  VECTOR6 lambda;
  lambda[0] = dTheta * (cosTheta - 1.0) / crossNorm;
  lambda[1] = dTheta * (cosTheta + 1.0) / crossNorm;
  const REAL radical = sqrt(alpha * alpha + 4 * dTheta * dTheta);
  lambda[2] = (alpha + radical) * 0.5;
  lambda[3] = (alpha - radical) * 0.5;
  lambda[4] = dTheta;
  lambda[5] = -dTheta;

  VECTOR6 q2, q3, q4, q5;
  const REAL denom0 = (fabs(dTheta + lambda[2]) > 0.0) ? dTheta + lambda[2] : 1.0;
  const REAL denom1 = (fabs(dTheta - lambda[2]) > 0.0) ? dTheta - lambda[2] : 1.0;
  const REAL denom2 = (fabs(dTheta + lambda[3]) > 0.0) ? dTheta + lambda[3] : 1.0;
  const REAL denom3 = (fabs(dTheta - lambda[3]) > 0.0) ? dTheta - lambda[3] : 1.0;
  q2 = (1.0 / denom0) * (q2theta + q4theta) + (1.0 / denom1) * (q3theta + q5theta);
  q3 = (1.0 / denom2) * (q2theta + q4theta) + (1.0 / denom3) * (q3theta + q5theta);
  q4 = 0.5 * (q3theta - q5theta);
  q5 = 0.5 * (q2theta - q4theta);

  q0.normalize(); q1.normalize(); q2.normalize();
  q3.normalize(); q4.normalize(); q5.normalize();

  MATRIX6 S;
  S.setZero();
  S(0,0) = S(1,1) = S(2,2) = 1.0 / e0.norm();
  S(3,3) = S(4,4) = S(5,5) = 1.0 / e1.norm();

  //S = [eye(3,3) / e0norm zeros(3,3); zeros(3,3) eye(3,3) / e1norm];
  q0 = S * q0; q1 = S * q1; q2 = S * q2;
  q3 = S * q3; q4 = S * q4; q5 = S * q5;

  q0.normalize(); q1.normalize(); q2.normalize();
  q3.normalize(); q4.normalize(); q5.normalize();

  MATRIX6 Q;
  Q.col(0) = q0; Q.col(1) = q1; Q.col(2) = q2;
  Q.col(3) = q3; Q.col(4) = q4; Q.col(5) = q5;
  
  const REAL scale = 0.5 * (1 / e0.dot(e0) + 1 / e1.dot(e1)); 
  lambda = lambda * scale;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();
  return H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_UNIT_BENDING::clampedHessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  return clampEigenvalues(hessian(E,theta0, inverted));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_UNIT_BENDING::clampedHessian(const MATRIX3x2& E, const REAL& theta0) const
{
  return clampEigenvalues(hessian(E,theta0));
  /*
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  MATRIX3x2 T;
  T.col(0) = t0;
  T.col(1) = t1;

  const REAL cosTheta = t0.dot(t1);
  const REAL crossNorm = (t0.cross(t1)).norm();

  if (fabs(crossNorm) < 1e-8 || fabs(acos(cosTheta)- theta0) < 1e-8)
  {
    // it will just be the outer product term, which is guaranteed semi-positive-definite
    VECTOR6 g = flatten(gradient(E));
    return _mu * g * g.transpose();
  }

  VECTOR6 lambda;
  lambda[0] = (cosTheta - 1.0) / crossNorm;
  lambda[1] = (cosTheta + 1.0) / crossNorm;

  const VECTOR3 crossed = t0.cross(t1).normalized();
  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);
  const VECTOR3 q2hat = t0 + e0cross;
  const VECTOR3 q3hat = t0 - e0cross;
  const VECTOR3 q4hat = t1 - e1cross;
  const VECTOR3 q5hat = t1 + e1cross;

  VECTOR6 q2theta;
  VECTOR6 q3theta;
  VECTOR6 q4theta;
  VECTOR6 q5theta;
  q2theta.setZero();
  q3theta.setZero();
  q4theta.setZero();
  q5theta.setZero();

  q2theta.block<3,1>(0,0) = q2hat;
  //q2theta[0] = q2hat[0]; q2theta[1] = q2hat[1]; q2theta[2] = q2hat[2];
  q3theta.block<3,1>(0,0) = q3hat;
  //q3theta[0] = q3hat[0]; q3theta[1] = q3hat[1]; q3theta[2] = q3hat[2];
  q4theta.block<3,1>(3,0) = q4hat;
  //q4theta[3] = q4hat[0]; q4theta[4] = q4hat[1]; q4theta[5] = q4hat[2];
  q5theta.block<3,1>(3,0) = q5hat;
  //q5theta[3] = q5hat[0]; q5theta[4] = q5hat[1]; q5theta[5] = q5hat[2];
  q2theta.normalize(); q3theta.normalize();
  q4theta.normalize(); q5theta.normalize();

  // segments calls bomb for some reason
  VECTOR6 q0;
  VECTOR6 q1;
  q0.block<3,1>(0,0) = crossed;
  //q0[0] = crossed[0]; q0[1] = crossed[1]; q0[2] = crossed[2];
  q0.block<3,1>(3,0) = crossed;
  //q0[3] = crossed[0]; q0[4] = crossed[1]; q0[5] = crossed[2];
  q1.block<3,1>(0,0) = crossed;
  //q1[0] = crossed[0]; q1[1] = crossed[1]; q1[2] = crossed[2];
  q1.block<3,1>(3,0) = -crossed;
  //q1[3] = -crossed[0]; q1[4] = -crossed[1]; q1[5] = -crossed[2];

  const VECTOR6 p = flatten(gradient(T));
  const REAL theta = angle(E);

  const REAL dTheta = _mu * (theta - theta0);
  const REAL ddTheta = _mu;

  // build scaled version from eigendecomposition
  const REAL alpha = ddTheta * p.dot(p);

  // MULT by a
  lambda[0] = dTheta * lambda[0];
  lambda[1] = dTheta * lambda[1];
  lambda[2] = (alpha + sqrt(alpha * alpha + 4 * dTheta * dTheta)) * 0.5;
  lambda[3] = (alpha - sqrt(alpha * alpha + 4 * dTheta * dTheta)) * 0.5;
  lambda[4] = dTheta;
  lambda[5] = -dTheta;

  VECTOR6 q2, q3, q4, q5;
  const REAL denom0 = (fabs(dTheta + lambda[2]) > 0.0) ? dTheta + lambda[2] : 1.0;
  const REAL denom1 = (fabs(dTheta - lambda[2]) > 0.0) ? dTheta - lambda[2] : 1.0;
  const REAL denom2 = (fabs(dTheta + lambda[3]) > 0.0) ? dTheta + lambda[3] : 1.0;
  const REAL denom3 = (fabs(dTheta - lambda[3]) > 0.0) ? dTheta - lambda[3] : 1.0;
  q2 = (1.0 / denom0) * (q2theta + q4theta) + (1.0 / denom1) * (q3theta + q5theta);
  q3 = (1.0 / denom2) * (q2theta + q4theta) + (1.0 / denom3) * (q3theta + q5theta);
  //q2 = (1 / (dTheta + lambda[2])) * (q2theta + q4theta) + (1 / (dTheta - lambda[2])) * (q3theta + q5theta);
  //q3 = (1 / (dTheta + lambda[3])) * (q2theta + q4theta) + (1 / (dTheta - lambda[3])) * (q3theta + q5theta);
  q4 = 0.5 * (q3theta - q5theta);
  q5 = 0.5 * (q2theta - q4theta);

  q0.normalize(); q1.normalize(); q2.normalize();
  q3.normalize(); q4.normalize(); q5.normalize();

  MATRIX6 S;
  S.setZero();
  S(0,0) = S(1,1) = S(2,2) = 1.0 / e0.norm();
  S(3,3) = S(4,4) = S(5,5) = 1.0 / e1.norm();

  //S = [eye(3,3) / e0norm zeros(3,3); zeros(3,3) eye(3,3) / e1norm];
  q0 = S * q0;
  q1 = S * q1;
  q2 = S * q2;
  q3 = S * q3;
  q4 = S * q4;
  q5 = S * q5;

  q0.normalize(); q1.normalize(); q2.normalize();
  q3.normalize(); q4.normalize(); q5.normalize();

  MATRIX6 Q;
  Q.col(0) = q0;
  Q.col(1) = q1;
  Q.col(2) = q2;
  Q.col(3) = q3;
  Q.col(4) = q4;
  Q.col(5) = q5;
  
  const REAL scale = 0.5 * (1 / e0.dot(e0) + 1 / e1.dot(e1)); 
  lambda = lambda * scale;
  for (unsigned int x = 0; x < 6; x++)
    lambda[x] = (lambda[x] > 0.0) ? lambda[x] : 0.0;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();
  
  return H;
  */
}

} // STRAND 
} // HOBAK
