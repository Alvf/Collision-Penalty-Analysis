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

#include "QUADRATIC_BENDING.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace STRAND {

QUADRATIC_BENDING::QUADRATIC_BENDING(const REAL& mu) :
  ISOTROPIC_BENDING(mu)
{
}

QUADRATIC_BENDING::QUADRATIC_BENDING(const REAL& mu, const REAL& theta0) :
  ISOTROPIC_BENDING(mu,theta0)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string QUADRATIC_BENDING::name() const
{ 
  return std::string("QUADRATIC_BENDING"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL QUADRATIC_BENDING::psi(const MATRIX3x2& E, const REAL& theta0) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const REAL diff = theta - theta0;

  return 0.5 * _mu * diff * diff;
  //return 0.5 * _mu * theta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL QUADRATIC_BENDING::psi(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const REAL sign = (inverted) ? 1.0 : -1.0;
  const REAL diff = theta - theta0;

  return 0.5 * _mu * diff * diff;
  //return 0.5 * _mu * theta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// gradient computation, once and for all
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 QUADRATIC_BENDING::gradient(const MATRIX3x2& E) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR2 u(1,0);
  const VECTOR2 v(0,1);

  const REAL I6 = (E * u).transpose() * (E * v);
  const REAL I5u = (E * u).transpose() * (E * u);
  const REAL I5v = (E * v).transpose() * (E * v);

  MATRIX2 anti, uut, vvt;
  anti << 0, 1, 1, 0;
  uut  << 1, 0, 0, 0;
  vvt  << 0, 0, 0, 1;

  const REAL alpha = I6 / (sqrt(I5u) * sqrt(I5v));
  const MATRIX3x2 P = (1 / (sqrt(I5u) * sqrt(I5v)) * E * anti -
                       I6 * (1 / (sqrt(I5v) * pow(I5u, 1.5)) * E * uut + 
                       1 / (sqrt(I5u) * pow(I5v,1.5)) * E * vvt));

  // this is slightly ugly, but seems to work
  const REAL alphaSq = alpha * alpha;
  const REAL denom = (alphaSq < 1.0) ? sqrt(1 - alpha * alpha) : 0.0;
  const REAL coeff = (fabs(denom) <= 0) ? -1.0 : -1.0 / denom;
  /*
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " alpha: " << alpha << endl;
  cout << " denom: " << denom << endl;
  cout << " coeff: " << coeff << endl;
  */
  return coeff * P;
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
MATRIX3x2 QUADRATIC_BENDING::PK1(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  //const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  //const REAL theta = acos(cosTheta);
  REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  REAL theta = acos(cosTheta);

  if (inverted)
  {
    theta = 2 * M_PI - theta;
    cosTheta = cos(theta);
  }

  const MATRIX3x2 Pfinal = _mu * (theta - theta0) * gradient(E);

  return Pfinal;
  //return _mu * (theta - theta0) * thetaPK1;
  //return _mu * thetaPK1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 QUADRATIC_BENDING::PK1(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  //const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  const REAL theta = acos(cosTheta);
  const MATRIX3x2 Pfinal = _mu * (theta - theta0) * gradient(E);
  /*
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " Pfinal: " << endl << Pfinal << endl;
  cout << " norm: " << Pfinal.squaredNorm() << endl;
  cout << " theta: " << theta << endl;
  cout << " E: " << endl<< E << endl;
  cout << " cosTheta: " << cosTheta << endl;
  cout << endl;
  */

  return Pfinal;
  //return _mu * (theta - theta0) * thetaPK1;
  //return _mu * thetaPK1;
}
#if 0
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR2 u(1,0);
  const VECTOR2 v(0,1);

  const REAL I6 = (E * u).transpose() * (E * v);
  const REAL I5u = (E * u).transpose() * (E * u);
  const REAL I5v = (E * v).transpose() * (E * v);

  MATRIX2 anti, uut, vvt;
  anti << 0, 1, 1, 0;
  uut  << 1, 0, 0, 0;
  vvt  << 0, 0, 0, 1;

  const REAL alpha = I6 / (sqrt(I5u) * sqrt(I5v));
  const MATRIX3x2 P = (1 / (sqrt(I5u) * sqrt(I5v)) * E * anti -
                       I6 * (1 / (sqrt(I5v) * pow(I5u, 1.5)) * E * uut + 
                       1 / (sqrt(I5u) * pow(I5v,1.5)) * E * vvt));

  /*
  const MATRIX3x2 thetaPK1 = (-1.0 / sqrt(1 - alpha * alpha)) * P;
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  return _mu * (theta - theta0) * thetaPK1;
  */

  // this is slightly ugly, but seems to work
  const REAL denom = sqrt(1 - alpha * alpha);
  const REAL coeff = (fabs(denom) <= 0) ? -1.0 : -1.0 / denom;
  const MATRIX3x2 thetaPK1 = coeff * P;
  
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const MATRIX3x2 Pfinal = _mu * (theta - theta0) * thetaPK1;
  cout << " Pfinal: " << Pfinal << endl;

  return Pfinal;
  //return _mu * (theta - theta0) * thetaPK1;
  //return _mu * thetaPK1;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_BENDING::hessian(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  const VECTOR3 crossed = (t0.cross(t1)).normalized();
  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);

  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();
  const REAL e0dot = e0.dot(e0);
  const REAL e1dot = e1.dot(e1);
  const REAL e0dotInv = 1.0 / e0dot;
  const REAL e1dotInv = 1.0 / e1dot;

  VECTOR6 lambda;
  lambda[2] = -e0dotInv;
  lambda[3] =  e0dotInv;
  lambda[4] = -e1dotInv;
  lambda[5] =  e1dotInv;

  const VECTOR3 q2hat = t0 + e0cross;
  const VECTOR3 q3hat = t0 - e0cross;
  const VECTOR3 q4hat = t1 - e1cross;
  const VECTOR3 q5hat = t1 + e1cross;

  VECTOR6 q2;
  VECTOR6 q3;
  VECTOR6 q4;
  VECTOR6 q5;
  q2.setZero();
  q3.setZero();
  q4.setZero();
  q5.setZero();

  //q2.segment(0,2) = q2hat;
  q2[0] = q2hat[0];
  q2[1] = q2hat[1];
  q2[2] = q2hat[2];
  //q3.segment(0,2) = q3hat;
  q3[0] = q3hat[0];
  q3[1] = q3hat[1];
  q3[2] = q3hat[2];
  //q4.segment(3,5) = q4hat;
  q4[3] = q4hat[0];
  q4[4] = q4hat[1];
  q4[5] = q4hat[2];
  //q5.segment(3,5) = q5hat;
  q5[3] = q5hat[0];
  q5[4] = q5hat[1];
  q5[5] = q5hat[2];

  const REAL cosTheta = e0.dot(e1) / (e0norm * e1norm);
  const REAL b = -(e0norm / e1norm - e1norm / e0norm) * cosTheta;
  const REAL rad = b*b + 4.0;
  const REAL ratio = (-b + sqrt(rad))/2;

  // segments calls bomb for some reason
  VECTOR6 q0;
  VECTOR6 q1;
  //q0.segment(0,2) = ratio * crossed;
  q0[0] = crossed[0] * ratio;
  q0[1] = crossed[1] * ratio;
  q0[2] = crossed[2] * ratio;
  //q0.segment(3,5) = crossed;
  q0[3] = crossed[0];
  q0[4] = crossed[1];
  q0[5] = crossed[2];
  //q1.segment(0,2) = crossed;
  q1[0] = crossed[0];
  q1[1] = crossed[1];
  q1[2] = crossed[2];
  //q1.segment(3,5) = -crossed * ratio;
  q1[3] = -crossed[0] * ratio;
  q1[4] = -crossed[1] * ratio;
  q1[5] = -crossed[2] * ratio;
 
  const REAL theta = acos(cosTheta);
  const REAL cot = 1.0 / tan(theta);
  const REAL csc = 1.0 / sin(theta);
  const REAL left = 0.5 * (e0dotInv + e1dotInv) * cot;
  const REAL right = sqrt(0.25 * (e0dotInv - e1dotInv) * (e0dotInv - e1dotInv) * cot * cot + 
                          csc * csc * e0dotInv * e1dotInv);
  lambda[0] = left - right;
  lambda[1] = left + right;

  q0.normalize();
  q1.normalize();
  q2.normalize();
  q3.normalize();
  q4.normalize();
  q5.normalize();

  MATRIX6 Q;
  Q.col(0) = q0;
  Q.col(1) = q1;
  Q.col(2) = q2;
  Q.col(3) = q3;
  Q.col(4) = q4;
  Q.col(5) = q5;

  const MATRIX6 thetaH = Q * lambda.asDiagonal() * Q.transpose();
  
  /*
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR2 u(1,0);
  const VECTOR2 v(0,1);

  const REAL I6 = (E * u).transpose() * (E * v);
  const REAL I5u = (E * u).transpose() * (E * u);
  const REAL I5v = (E * v).transpose() * (E * v);

  MATRIX2 anti, uut, vvt;
  anti << 0, 1, 1, 0;
  uut  << 1, 0, 0, 0;
  vvt  << 0, 0, 0, 1;

  const REAL alpha = I6 / (sqrt(I5u) * sqrt(I5v));
  const MATRIX3x2 P = (1 / (sqrt(I5u) * sqrt(I5v)) * E * anti -
                       I6 * (1 / (sqrt(I5v) * pow(I5u, 1.5)) * E * uut + 
                       1 / (sqrt(I5u) * pow(I5v,1.5)) * E * vvt));
  const VECTOR6 thetaP = (-1.0 / sqrt(1 - alpha * alpha)) * flatten(P);
  */
  const VECTOR6 thetaP = flatten(gradient(E));

  return _mu * (thetaP * thetaP.transpose()) + _mu * (theta - theta0) * thetaH;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_BENDING::hessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  const VECTOR3 crossed = (t0.cross(t1)).normalized();
  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);

  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();
  const REAL e0dot = e0.dot(e0);
  const REAL e1dot = e1.dot(e1);
  const REAL e0dotInv = 1.0 / e0dot;
  const REAL e1dotInv = 1.0 / e1dot;

  VECTOR6 lambda;
  lambda[2] = -e0dotInv;
  lambda[3] =  e0dotInv;
  lambda[4] = -e1dotInv;
  lambda[5] =  e1dotInv;

  const VECTOR3 q2hat = t0 + e0cross;
  const VECTOR3 q3hat = t0 - e0cross;
  const VECTOR3 q4hat = t1 - e1cross;
  const VECTOR3 q5hat = t1 + e1cross;

  VECTOR6 q2;
  VECTOR6 q3;
  VECTOR6 q4;
  VECTOR6 q5;
  q2.setZero();
  q3.setZero();
  q4.setZero();
  q5.setZero();

  //q2.segment(0,2) = q2hat;
  q2[0] = q2hat[0];
  q2[1] = q2hat[1];
  q2[2] = q2hat[2];
  //q3.segment(0,2) = q3hat;
  q3[0] = q3hat[0];
  q3[1] = q3hat[1];
  q3[2] = q3hat[2];
  //q4.segment(3,5) = q4hat;
  q4[3] = q4hat[0];
  q4[4] = q4hat[1];
  q4[5] = q4hat[2];
  //q5.segment(3,5) = q5hat;
  q5[3] = q5hat[0];
  q5[4] = q5hat[1];
  q5[5] = q5hat[2];

  const REAL cosTheta = e0.dot(e1) / (e0norm * e1norm);
  const REAL b = -(e0norm / e1norm - e1norm / e0norm) * cosTheta;
  const REAL rad = b*b + 4.0;
  const REAL ratio = (-b + sqrt(rad))/2;

  // segments calls bomb for some reason
  VECTOR6 q0;
  VECTOR6 q1;
  //q0.segment(0,2) = ratio * crossed;
  q0[0] = crossed[0] * ratio;
  q0[1] = crossed[1] * ratio;
  q0[2] = crossed[2] * ratio;
  //q0.segment(3,5) = crossed;
  q0[3] = crossed[0];
  q0[4] = crossed[1];
  q0[5] = crossed[2];
  //q1.segment(0,2) = crossed;
  q1[0] = crossed[0];
  q1[1] = crossed[1];
  q1[2] = crossed[2];
  //q1.segment(3,5) = -crossed * ratio;
  q1[3] = -crossed[0] * ratio;
  q1[4] = -crossed[1] * ratio;
  q1[5] = -crossed[2] * ratio;
 
  const REAL theta = acos(cosTheta);
  const REAL cot = 1.0 / tan(theta);
  const REAL csc = 1.0 / sin(theta);
  const REAL left = 0.5 * (e0dotInv + e1dotInv) * cot;
  const REAL right = sqrt(0.25 * (e0dotInv - e1dotInv) * (e0dotInv - e1dotInv) * cot * cot + 
                          csc * csc * e0dotInv * e1dotInv);
  lambda[0] = left - right;
  lambda[1] = left + right;

  q0.normalize();
  q1.normalize();
  q2.normalize();
  q3.normalize();
  q4.normalize();
  q5.normalize();

  MATRIX6 Q;
  Q.col(0) = q0;
  Q.col(1) = q1;
  Q.col(2) = q2;
  Q.col(3) = q3;
  Q.col(4) = q4;
  Q.col(5) = q5;

  const MATRIX6 thetaH = Q * lambda.asDiagonal() * Q.transpose();
  
  /*
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR2 u(1,0);
  const VECTOR2 v(0,1);

  const REAL I6 = (E * u).transpose() * (E * v);
  const REAL I5u = (E * u).transpose() * (E * u);
  const REAL I5v = (E * v).transpose() * (E * v);

  MATRIX2 anti, uut, vvt;
  anti << 0, 1, 1, 0;
  uut  << 1, 0, 0, 0;
  vvt  << 0, 0, 0, 1;

  const REAL alpha = I6 / (sqrt(I5u) * sqrt(I5v));
  const MATRIX3x2 P = (1 / (sqrt(I5u) * sqrt(I5v)) * E * anti -
                       I6 * (1 / (sqrt(I5v) * pow(I5u, 1.5)) * E * uut + 
                       1 / (sqrt(I5u) * pow(I5v,1.5)) * E * vvt));
  const VECTOR6 thetaP = (-1.0 / sqrt(1 - alpha * alpha)) * flatten(P);
  */
  const VECTOR6 thetaP = flatten(gradient(E));
  const REAL sign = (inverted) ? 1.0 : -1.0;

  return _mu * (thetaP * thetaP.transpose()) + _mu * (theta - theta0) * thetaH;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_BENDING::clampedHessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  const VECTOR3 crossed = (t0.cross(t1)).normalized();

  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);

  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();
  const REAL cosTheta = e0.dot(e1) / (e0norm * e1norm);
  const REAL theta = acos(cosTheta);


  const REAL sign = (inverted) ? 1.0 : -1.0;
  if (fabs(theta - theta0) < 1e-8)
  {
    // it will just be the outer product term, which is guaranteed semi-positive-definite
    VECTOR6 g = flatten(gradient(E));
    return _mu * g * g.transpose();
  }

  const REAL e0dot = e0.dot(e0);
  const REAL e1dot = e1.dot(e1);
  const REAL e0dotInv = 1.0 / e0dot;
  const REAL e1dotInv = 1.0 / e1dot;

  VECTOR6 lambdaTheta;
  lambdaTheta[2] = -e0dotInv;
  lambdaTheta[3] =  e0dotInv;
  lambdaTheta[4] = -e1dotInv;
  lambdaTheta[5] =  e1dotInv;

  const VECTOR3 q2hat = t0 + e0cross;
  const VECTOR3 q3hat = t0 - e0cross;
  const VECTOR3 q4hat = t1 - e1cross;
  const VECTOR3 q5hat = t1 + e1cross;

  VECTOR6 qTheta2;
  VECTOR6 qTheta3;
  VECTOR6 qTheta4;
  VECTOR6 qTheta5;
  qTheta2.setZero();
  qTheta3.setZero();
  qTheta4.setZero();
  qTheta5.setZero();

  //q2.segment(0,2) = q2hat;
  qTheta2[0] = q2hat[0];
  qTheta2[1] = q2hat[1];
  qTheta2[2] = q2hat[2];
  //q3.segment(0,2) = q3hat;
  qTheta3[0] = q3hat[0];
  qTheta3[1] = q3hat[1];
  qTheta3[2] = q3hat[2];
  //q4.segment(3,5) = q4hat;
  qTheta4[3] = q4hat[0];
  qTheta4[4] = q4hat[1];
  qTheta4[5] = q4hat[2];
  //q5.segment(3,5) = q5hat;
  qTheta5[3] = q5hat[0];
  qTheta5[4] = q5hat[1];
  qTheta5[5] = q5hat[2];

  const REAL d = -(e0norm / e1norm - e1norm / e0norm) * cosTheta;
  const REAL rad = d*d + 4.0;
  const REAL ratio = (-d + sqrt(rad))/2;

  // segments calls bomb for some reason
  VECTOR6 qTheta0;
  VECTOR6 qTheta1;
  //q0.segment(0,2) = ratio * crossed;
  qTheta0[0] = crossed[0] * ratio;
  qTheta0[1] = crossed[1] * ratio;
  qTheta0[2] = crossed[2] * ratio;
  //q0.segment(3,5) = crossed;
  qTheta0[3] = crossed[0];
  qTheta0[4] = crossed[1];
  qTheta0[5] = crossed[2];
  //q1.segment(0,2) = crossed;
  qTheta1[0] = crossed[0];
  qTheta1[1] = crossed[1];
  qTheta1[2] = crossed[2];
  //q1.segment(3,5) = -crossed * ratio;
  qTheta1[3] = -crossed[0] * ratio;
  qTheta1[4] = -crossed[1] * ratio;
  qTheta1[5] = -crossed[2] * ratio;
 
  const REAL cot = 1.0 / tan(theta);
  const REAL csc = 1.0 / sin(theta);
  const REAL left = 0.5 * (e0dotInv + e1dotInv) * cot;
  const REAL right = sqrt(0.25 * (e0dotInv - e1dotInv) * (e0dotInv - e1dotInv) * cot * cot + 
                          csc * csc * e0dotInv * e1dotInv);
  lambdaTheta[0] = left - right;
  lambdaTheta[1] = left + right;

  qTheta0.normalize();
  qTheta1.normalize();
  qTheta2.normalize();
  qTheta3.normalize();
  qTheta4.normalize();
  qTheta5.normalize();

  const REAL c = fabs(theta - theta0) > 0 ? 1.0 / (theta - theta0) : 1.0;
  const REAL a = e0dot;
  const REAL b = e1dot;
  const REAL a2 = a * a;
  const REAL b2 = b * b;
  const REAL c2 = c * c;

  //alpha = sqrt((a^2 + b^2) * (4 + c^2) + a * b * (2 * c^2 - 8));
  const REAL alpha = sqrt((a2 + b2) * (4 + c2) + a * b * (2.0 * c2 - 8.0));
  //eta = (c^3 * (a + b)^3 + 4 * c * (a - b)^2 * (a + b));
  const REAL eta = (c * c2 * pow((a + b),3) + 4.0 * c * (a - b) * (a - b) * (a + b));
  //left = 0.5 * c^2 * (a + b)^2 + (a^2 + b^2) + 2 * a * b;
  const REAL l = 0.5 * c2 * (a + b) * (a + b) + (a2 + b2) + 2 * a * b;
  //right = eta/(2* alpha);
  const REAL r = eta / (2.0 * alpha);

  VECTOR6 lambda;
  lambda[0] = lambdaTheta[0];
  lambda[1] = lambdaTheta[1];
  const VECTOR6 q0 = qTheta0;
  const VECTOR6 q1 = qTheta1;

  const REAL gamma0 = l + r;
  lambda[2] = (1.0 / (4.0 * a * b)) * ((alpha + c * (a + b)) - 2.0 * sqrt(gamma0));
  lambda[3] = (1.0 / (4.0 * a * b)) * ((alpha + c * (a + b)) + 2.0 * sqrt(gamma0));

  const REAL gamma1 = l - r;
  lambda[4] = (1.0 / (4.0 * a * b)) * ((-alpha + c * (a + b)) - 2.0 * sqrt(gamma1));
  lambda[5] = (1.0 / (4.0 * a * b)) * ((-alpha + c * (a + b)) + 2.0 * sqrt(gamma1));
  
  const bool verbose = false;
  if (verbose)
  {
    cout << " 4ab: " << 4.0 * a * b << endl;
    cout << " gamma: " << gamma0 << " " << gamma1 << endl;
    cout << " alpha: " << alpha << endl;
    cout << " a: " << a << endl;
    cout << " b: " << b << endl;
    cout << " c: " << c << endl;
  }

  //Q = zeros(6,6);
  MATRIX6 Q;
  Q.col(0) = qTheta0;
  Q.col(1) = qTheta1;
  
  //V = [lambda0 lambda1 lambda2 lambda3 lambda4 lambda5];

  //for i = 3:6
  for (unsigned int i = 2; i < 6; i++)
  {
    const REAL l = lambda[i];
    const REAL denom2 = -(e0dot * l + 1);
    const REAL denom3 = (e0dot * l - 1);
    const REAL denom4 = -(e1dot *  l + 1);
    const REAL denom5 = (e1dot * l - 1);
#if 0
    const REAL coeff2 = (fabs(denom2) > 0.0) ? -e0norm / denom2 : 0.0;
    const REAL coeff3 = (fabs(denom3) > 0.0) ? e0norm / denom3 : 0.0;
    const REAL coeff4 = (fabs(denom4) > 0.0) ? -e1norm / denom4 : 0.0;
    const REAL coeff5 = (fabs(denom5) > 0.0) ? e1norm / denom5 : 0.0;
#else
    REAL coeff2 = -e0norm / (e0dot * l + 1);
    REAL coeff3 = e0norm / (e0dot * l - 1);
    REAL coeff4 = -e1norm / (e1dot *  l + 1);
    REAL coeff5 = e1norm / (e1dot * l - 1);

    if (fabs(denom2) < 1e-8)
    {
      coeff2 = -coeff3;
      coeff3 = 0;
    }
    if (fabs(denom3) < 1e-8)
    {
      coeff3 = -coeff2;
      coeff2 = 0;
    }
    if (fabs(denom4) < 1e-8)
    {
      coeff4 = coeff5;
      coeff5 = 0;
    }
    if (fabs(denom5) < 1e-8)
    {
      coeff5 = coeff4;
      coeff4 = 0;
    }
#endif
    const VECTOR6 q = coeff2 * qTheta2 + coeff3 * qTheta3 + coeff4 * qTheta4 + coeff5 * qTheta5;

    if (verbose)
    {
      cout << " coeffs: " << coeff2 << " " << coeff3 << " " << coeff4 << " " << coeff5 << endl;
      cout << " denoms: " << denom2 << " " << denom3 << " " << denom4 << " " << denom5 << endl;
      cout << " lambda: " << l << endl;
      cout << " dots: " << e0dot << " " << e1dot << endl;
    }
    Q.col(i) = q.normalized();
  }
  if (verbose)
  {
    cout << " lambda: " << lambda << endl;
    cout << " Q: " << endl << Q << endl;
    cout << "theta - theta0: " << theta - theta0 << endl; 
  }

  lambda *= _mu * (theta - theta0);
  for (int i = 0; i < 6; i++)
    lambda[i] = (lambda[i] > 0.0) ? lambda[i] : 0.0;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();

  //cout << " H: " << endl << H << endl;

  return H;
  //return Q * lambda.asDiagonal() * Q.transpose();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_BENDING::clampedHessian(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const VECTOR3 t0 = e0 / e0.norm();
  const VECTOR3 t1 = e1 / e1.norm();

  const VECTOR3 crossed = (t0.cross(t1)).normalized();

  const VECTOR3 e0cross = t0.cross(crossed);
  const VECTOR3 e1cross = t1.cross(crossed);

  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();
  const REAL cosTheta = e0.dot(e1) / (e0norm * e1norm);
  const REAL theta = acos(cosTheta);

  if (fabs(theta - theta0) < 1e-8)
  {
    // it will just be the outer product term, which is guaranteed semi-positive-definite
    VECTOR6 g = flatten(gradient(E));
    return _mu * g * g.transpose();
  }

  const REAL e0dot = e0.dot(e0);
  const REAL e1dot = e1.dot(e1);
  const REAL e0dotInv = 1.0 / e0dot;
  const REAL e1dotInv = 1.0 / e1dot;

  VECTOR6 lambdaTheta;
  lambdaTheta[2] = -e0dotInv;
  lambdaTheta[3] =  e0dotInv;
  lambdaTheta[4] = -e1dotInv;
  lambdaTheta[5] =  e1dotInv;

  const VECTOR3 q2hat = t0 + e0cross;
  const VECTOR3 q3hat = t0 - e0cross;
  const VECTOR3 q4hat = t1 - e1cross;
  const VECTOR3 q5hat = t1 + e1cross;

  VECTOR6 qTheta2;
  VECTOR6 qTheta3;
  VECTOR6 qTheta4;
  VECTOR6 qTheta5;
  qTheta2.setZero();
  qTheta3.setZero();
  qTheta4.setZero();
  qTheta5.setZero();

  //q2.segment(0,2) = q2hat;
  qTheta2[0] = q2hat[0];
  qTheta2[1] = q2hat[1];
  qTheta2[2] = q2hat[2];
  //q3.segment(0,2) = q3hat;
  qTheta3[0] = q3hat[0];
  qTheta3[1] = q3hat[1];
  qTheta3[2] = q3hat[2];
  //q4.segment(3,5) = q4hat;
  qTheta4[3] = q4hat[0];
  qTheta4[4] = q4hat[1];
  qTheta4[5] = q4hat[2];
  //q5.segment(3,5) = q5hat;
  qTheta5[3] = q5hat[0];
  qTheta5[4] = q5hat[1];
  qTheta5[5] = q5hat[2];

  const REAL d = -(e0norm / e1norm - e1norm / e0norm) * cosTheta;
  const REAL rad = d*d + 4.0;
  const REAL ratio = (-d + sqrt(rad))/2;

  // segments calls bomb for some reason
  VECTOR6 qTheta0;
  VECTOR6 qTheta1;
  //q0.segment(0,2) = ratio * crossed;
  qTheta0[0] = crossed[0] * ratio;
  qTheta0[1] = crossed[1] * ratio;
  qTheta0[2] = crossed[2] * ratio;
  //q0.segment(3,5) = crossed;
  qTheta0[3] = crossed[0];
  qTheta0[4] = crossed[1];
  qTheta0[5] = crossed[2];
  //q1.segment(0,2) = crossed;
  qTheta1[0] = crossed[0];
  qTheta1[1] = crossed[1];
  qTheta1[2] = crossed[2];
  //q1.segment(3,5) = -crossed * ratio;
  qTheta1[3] = -crossed[0] * ratio;
  qTheta1[4] = -crossed[1] * ratio;
  qTheta1[5] = -crossed[2] * ratio;
 
  const REAL cot = 1.0 / tan(theta);
  const REAL csc = 1.0 / sin(theta);
  const REAL left = 0.5 * (e0dotInv + e1dotInv) * cot;
  const REAL right = sqrt(0.25 * (e0dotInv - e1dotInv) * (e0dotInv - e1dotInv) * cot * cot + 
                          csc * csc * e0dotInv * e1dotInv);
  lambdaTheta[0] = left - right;
  lambdaTheta[1] = left + right;

  qTheta0.normalize();
  qTheta1.normalize();
  qTheta2.normalize();
  qTheta3.normalize();
  qTheta4.normalize();
  qTheta5.normalize();

  const REAL c = fabs(theta - theta0) > 0 ? 1.0 / (theta - theta0) : 1.0;
  const REAL a = e0dot;
  const REAL b = e1dot;
  const REAL a2 = a * a;
  const REAL b2 = b * b;
  const REAL c2 = c * c;

  //alpha = sqrt((a^2 + b^2) * (4 + c^2) + a * b * (2 * c^2 - 8));
  const REAL alpha = sqrt((a2 + b2) * (4 + c2) + a * b * (2.0 * c2 - 8.0));
  //eta = (c^3 * (a + b)^3 + 4 * c * (a - b)^2 * (a + b));
  const REAL eta = (c * c2 * pow((a + b),3) + 4.0 * c * (a - b) * (a - b) * (a + b));
  //left = 0.5 * c^2 * (a + b)^2 + (a^2 + b^2) + 2 * a * b;
  const REAL l = 0.5 * c2 * (a + b) * (a + b) + (a2 + b2) + 2 * a * b;
  //right = eta/(2* alpha);
  const REAL r = eta / (2.0 * alpha);

  VECTOR6 lambda;
  lambda[0] = lambdaTheta[0];
  lambda[1] = lambdaTheta[1];
  const VECTOR6 q0 = qTheta0;
  const VECTOR6 q1 = qTheta1;

  const REAL gamma0 = l + r;
  lambda[2] = (1.0 / (4.0 * a * b)) * ((alpha + c * (a + b)) - 2.0 * sqrt(gamma0));
  lambda[3] = (1.0 / (4.0 * a * b)) * ((alpha + c * (a + b)) + 2.0 * sqrt(gamma0));

  const REAL gamma1 = l - r;
  lambda[4] = (1.0 / (4.0 * a * b)) * ((-alpha + c * (a + b)) - 2.0 * sqrt(gamma1));
  lambda[5] = (1.0 / (4.0 * a * b)) * ((-alpha + c * (a + b)) + 2.0 * sqrt(gamma1));
  
  const bool verbose = false;
  if (verbose)
  {
    cout << " 4ab: " << 4.0 * a * b << endl;
    cout << " gamma: " << gamma0 << " " << gamma1 << endl;
    cout << " alpha: " << alpha << endl;
    cout << " a: " << a << endl;
    cout << " b: " << b << endl;
    cout << " c: " << c << endl;
  }

  //Q = zeros(6,6);
  MATRIX6 Q;
  Q.col(0) = qTheta0;
  Q.col(1) = qTheta1;
  
  //V = [lambda0 lambda1 lambda2 lambda3 lambda4 lambda5];

  //for i = 3:6
  for (unsigned int i = 2; i < 6; i++)
  {
    const REAL l = lambda[i];
    const REAL denom2 = -(e0dot * l + 1);
    const REAL denom3 = (e0dot * l - 1);
    const REAL denom4 = -(e1dot *  l + 1);
    const REAL denom5 = (e1dot * l - 1);
#if 0
    const REAL coeff2 = (fabs(denom2) > 0.0) ? -e0norm / denom2 : 0.0;
    const REAL coeff3 = (fabs(denom3) > 0.0) ? e0norm / denom3 : 0.0;
    const REAL coeff4 = (fabs(denom4) > 0.0) ? -e1norm / denom4 : 0.0;
    const REAL coeff5 = (fabs(denom5) > 0.0) ? e1norm / denom5 : 0.0;
#else
    REAL coeff2 = -e0norm / (e0dot * l + 1);
    REAL coeff3 = e0norm / (e0dot * l - 1);
    REAL coeff4 = -e1norm / (e1dot *  l + 1);
    REAL coeff5 = e1norm / (e1dot * l - 1);

    if (fabs(denom2) < 1e-8)
    {
      coeff2 = -coeff3;
      coeff3 = 0;
    }
    if (fabs(denom3) < 1e-8)
    {
      coeff3 = -coeff2;
      coeff2 = 0;
    }
    if (fabs(denom4) < 1e-8)
    {
      coeff4 = coeff5;
      coeff5 = 0;
    }
    if (fabs(denom5) < 1e-8)
    {
      coeff5 = coeff4;
      coeff4 = 0;
    }
#endif
    const VECTOR6 q = coeff2 * qTheta2 + coeff3 * qTheta3 + coeff4 * qTheta4 + coeff5 * qTheta5;

    if (verbose)
    {
      cout << " coeffs: " << coeff2 << " " << coeff3 << " " << coeff4 << " " << coeff5 << endl;
      cout << " denoms: " << denom2 << " " << denom3 << " " << denom4 << " " << denom5 << endl;
      cout << " lambda: " << l << endl;
      cout << " dots: " << e0dot << " " << e1dot << endl;
    }
    Q.col(i) = q.normalized();
  }
  if (verbose)
  {
    cout << " lambda: " << lambda << endl;
    cout << " Q: " << endl << Q << endl;
    cout << "theta - theta0: " << theta - theta0 << endl; 
  }

  lambda *= _mu * (theta - theta0);
  for (int i = 0; i < 6; i++)
    lambda[i] = (lambda[i] > 0.0) ? lambda[i] : 0.0;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();

  //cout << " H: " << endl << H << endl;

  return H;
  //return Q * lambda.asDiagonal() * Q.transpose();
}

} // STRAND 
} // HOBAK
