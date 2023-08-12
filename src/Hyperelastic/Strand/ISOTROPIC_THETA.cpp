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

#include "ISOTROPIC_THETA.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace STRAND {

ISOTROPIC_THETA::ISOTROPIC_THETA(const REAL& mu) :
  ISOTROPIC_BENDING(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string ISOTROPIC_THETA::name() const
{ 
  return std::string("ISOTROPIC_THETA"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL ISOTROPIC_THETA::psi(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  return _mu * acos(cosTheta);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL ISOTROPIC_THETA::psi(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  REAL theta = acos(cosTheta);
  if (inverted)
    theta = 2.0 * M_PI - theta;
  return _mu * theta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 ISOTROPIC_THETA::PK1(const MATRIX3x2& E, const REAL& theta0) const
{
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
  MATRIX3x2 P = (1 / (sqrt(I5u) * sqrt(I5v)) * E * anti -
                 I6 * (1 / (sqrt(I5v) * pow(I5u, 1.5)) * E * uut + 
                       1 / (sqrt(I5u) * pow(I5v,1.5)) * E * vvt));

  /*
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " P: " << endl << P << endl;
  cout << " theta: " << acos(I6 / (sqrt(I5u) * sqrt(I5v))) << endl;
  cout << " alpha: " << alpha << endl;
  */

  return _mu * (-1.0 / sqrt(1 - alpha * alpha)) * P;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 ISOTROPIC_THETA::PK1(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " NOT IMPLEMENTED. " << endl;
  exit(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 ISOTROPIC_THETA::hessian(const MATRIX3x2& E, const REAL& theta0) const
{
  //u = [1 0]';
  //v = [0 1]';
  //I6 = (F * u)' * (F * v);
  //I5u = (F * u)' * (F * u);
  //I5v = (F * v)' * (F * v);
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL I6 = e0.dot(e1);
  const REAL I5u = e0.dot(e0);
  const REAL I5v = e1.dot(e1);
 
  MATRIX2 anti, uut, vvt;
  anti << 0, 1, 1, 0;
  uut  << 1, 0, 0, 0;
  vvt  << 0, 0, 0, 1;

  MATRIX6 FF, FFu, FFv;
  //FF = [zeros(3,3) eye(3,3);
  //      eye(3,3) zeros(3,3)];
  FF.setZero();
  FF.block<3,3>(3,0) = MATRIX3::Identity();
  FF.block<3,3>(0,3) = MATRIX3::Identity();
  
  //FFu = [eye(3,3) zeros(3,3);
  //       zeros(3,3) zeros(3,3)];
  FFu.setZero();
  FFu.block<3,3>(0,0) = MATRIX3::Identity();

  //FFv = [zeros(3,3) zeros(3,3);
  //       zeros(3,3) eye(3,3)];
  FFv.setZero();
  FFv.block<3,3>(3,3) = MATRIX3::Identity();

  const MATRIX3x2 Fu = E * uut;
  const MATRIX3x2 Fv = E * vvt;
  const MATRIX3x2 Fa = E * anti;
  const VECTOR6 fu = flatten(Fu);
  const VECTOR6 fv = flatten(Fv);
  const VECTOR6 fa = flatten(Fa);
  
  const REAL alpha = I6 / (sqrt(I5u) * sqrt(I5v));
  const MATRIX3x2 P = (1 / (sqrt(I5u) * sqrt(I5v)) * E * anti -
                      I6 * (1 / (sqrt(I5v) * pow(I5u, 1.5)) * E * uut + 
                            1 / (sqrt(I5u) * pow(I5v,1.5)) * E * vvt));
  const VECTOR6 p = flatten(P);

  //beast = - 1 / (I5u^1.5 * sqrt(I5v)) * (fa * fu' + fu * fa') ...
  //    - 1 / (sqrt(I5u) * I5v^1.5) * (fa * fv' + fv * fa') ...
  //    + I6 / (I5u^1.5 * I5v^1.5) * (fu * fv' + fv * fu') ...
  //    + (3 * I6 / (I5u^2.5 * sqrt(I5v))) * fu * fu' ...
  //    + (3 * I6 / (sqrt(I5u) * I5v^2.5)) * fv * fv' ...
  //    + 1 / (sqrt(I5u) * sqrt(I5v)) * FF ...
  //    - (I6 / (I5u^1.5 * sqrt(I5v))) * FFu ...
  //    - (I6 / (sqrt(I5u) * I5v^1.5)) * FFv;
  const REAL I5u15 = pow(I5u, 1.5);
  const REAL I5u25 = pow(I5u, 2.5);
  const REAL I5v15 = pow(I5v, 1.5);
  const REAL I5v25 = pow(I5v, 2.5);
  const MATRIX6 beast = - 1.0 / (I5u15 * sqrt(I5v)) * (fa * fu.transpose() + fu * fa.transpose())
                        - 1.0 / (sqrt(I5u) * I5v15) * (fa * fv.transpose() + fv * fa.transpose())
                        + I6 / (I5u15 * I5v15) * (fu * fv.transpose() + fv * fu.transpose())
                        + (3.0 * I6 / (I5u25 * sqrt(I5v))) * fu * fu.transpose()
                        + (3.0 * I6 / (sqrt(I5u) * I5v25)) * fv * fv.transpose()
                        + 1.0 / (sqrt(I5u) * sqrt(I5v)) * FF
                        - (I6 / (I5u15 * sqrt(I5v))) * FFu
                        - (I6 / (sqrt(I5u) * I5v15)) * FFv;
  //H = -(1 / sqrt(1 - alpha * alpha)) * beast - (alpha / (1 - alpha * alpha)^1.5) * p * p';
  const MATRIX6 H = -(1.0 / sqrt(1 - alpha * alpha)) * beast - (alpha / pow((1 - alpha * alpha),1.5)) * p * p.transpose();

  return _mu * H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 ISOTROPIC_THETA::hessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " NOT IMPLEMENTED. " << endl;
  exit(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 ISOTROPIC_THETA::clampedHessian(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);

  //e0normed = e0 / norm(e0);
  const VECTOR3 t0 = e0 / e0.norm();
  //e1normed = e1 / norm(e1);
  const VECTOR3 t1 = e1 / e1.norm();

  //crossed = cross(e0normed, e1normed);
  const VECTOR3 crossed = (t0.cross(t1)).normalized();
  //crossed = crossed / norm(crossed);
  //crossed = crossed / crossed.norm();

  //e0cross = cross(e0normed, crossed);
  const VECTOR3 e0cross = t0.cross(crossed);
  //e1cross = cross(e1normed, crossed);
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

  // do the filtering
  for (int i = 0; i < 6; i++)
    lambda[i] = (lambda[i] > 0.0) ? lambda[i] : 0.0;

  return _mu * (Q * lambda.asDiagonal() * Q.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 ISOTROPIC_THETA::clampedHessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " NOT IMPLEMENTED. " << endl;
  exit(0);
}

} // STRAND 
} // HOBAK
