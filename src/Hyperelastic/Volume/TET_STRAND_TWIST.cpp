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
#include "TET_STRAND_TWIST.h"
#include "MATRIX_UTIL.h"

#include <iostream>
#include <cmath>
#include "ext/solvePoly/poly34.h"

using namespace std;

namespace HOBAK {
namespace VOLUME {

using namespace std;

TET_STRAND_TWIST::TET_STRAND_TWIST(const REAL& mu, const REAL& theta0) :
    _mu(mu), _theta0(theta0)
{
  //cout << " theta0: " << _theta0 << endl;
  _zeroGuard = 1e-8;
}

std::string TET_STRAND_TWIST::name() const
{
  return "TET_STRAND_TWIST";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL TET_STRAND_TWIST::psi(const MATRIX3& F) const
{
  VECTOR3 u(1,0,0);
  VECTOR3 v(0,1,0);
  REAL i6 = (F*u).dot(F*v);
  REAL i5u = invariant5(F,u);
  REAL i5v = invariant5(F,v);

  REAL i3 = invariant3(F);
  REAL sgn = i3 > 0 ? -1.0 : 1.0;

  REAL diff = sgn*acos(i6/sqrt(i5u*i5v)) - _theta0;

  const REAL angle = acos(i6/sqrt(i5u*i5v));
  //cout << " theta0: " << _theta0 << " angle: " << angle << " diff: " << diff << endl;
  return _mu*diff*diff;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 TET_STRAND_TWIST::PK1(const MATRIX3& F) const
{
  VECTOR3 u(1,0,0);
  VECTOR3 v(0,1,0);
  REAL i6 = (F*u).dot(F*v);
  REAL i5u = invariant5(F,u);
  REAL i5v = invariant5(F,v);

  MATRIX3 Pi6PF = F*(u*v.transpose() + v*u.transpose());
  MATRIX3 Pi5uPF = 2*F*u*u.transpose();
  MATRIX3 Pi5vPF = 2*F*v*v.transpose();

  REAL i3 = invariant3(F);
  if(i3<_zeroGuard && i3 > -_zeroGuard) return MATRIX3::Zero(); //return 0 if NaNs are approaching.
  REAL sgn = i3 > 0 ? -1.0 : 1.0;

  REAL sqrtuv = sqrt(i5u*i5v);
  //cout << " sqrtuv: " << sqrtuv << endl;
  if(i6*i6/(i5u*i5v) - 1 < _zeroGuard && i6*i6/(i5u*i5v) - 1 > -_zeroGuard) return MATRIX3::Zero(); //return 0 if NaNs approaching 
  REAL sqrt1Diff = sqrt(1 - i6*i6/(i5u*i5v));
  //cout << " sqrt1Diff: " << sqrt1Diff<< endl;
  REAL diff = sgn*acos(i6/sqrtuv) - _theta0;
  /*
  cout << " i6: " << i6 << endl;
  cout << " acos: " << acos(i6 / sqrtuv) << endl;
  cout << " sgn: " << sgn << endl;
  cout << " theta0: " << _theta0 << endl;
  cout << " diff: " << diff << endl;
  */

  REAL DPsiDi6 = -2*_mu*sgn*diff/(sqrtuv*sqrt1Diff);
  REAL DPsiDi5u = i5v*i6*_mu*sgn*diff/(sqrtuv*sqrtuv*sqrtuv*sqrt1Diff);
  REAL DPsiDi5v = i5u*i6*_mu*sgn*diff/(sqrtuv*sqrtuv*sqrtuv*sqrt1Diff);

  /*
  cout << "twist F: " << endl << F << endl;
  MATRIX3 result = DPsiDi6*Pi6PF + DPsiDi5u*Pi5uPF + DPsiDi5v*Pi5vPF;
  cout << "twist PK1: " << endl << result << endl;
  MATRIX3 term1 = DPsiDi6*Pi6PF;
  MATRIX3 term2 = DPsiDi5u*Pi5uPF;
  MATRIX3 term3 = DPsiDi5v*Pi5vPF;
  cout << "term 1: " << endl << term1 << endl;
  cout << "term 2: " << endl << term2 << endl;
  cout << "term 3: " << endl << term3 << endl;
  cout << "DPsiDi6: " << DPsiDi6 << endl;
  cout << "DPsiDi5u: "<< DPsiDi5u << endl;
  cout << "DPsiDi5v: "<< DPsiDi5v << endl;
  */

  return DPsiDi6*Pi6PF + DPsiDi5u*Pi5uPF + DPsiDi5v*Pi5vPF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 TET_STRAND_TWIST::hessian(const MATRIX3& F) const
{
  VECTOR3 u(1,0,0);
  VECTOR3 v(0,1,0);
  REAL i6 = (F*u).dot(F*v);
  REAL i5u = invariant5(F,u);
  REAL i5v = invariant5(F,v);

  MATRIX3 Pi6PF = F*(u*v.transpose() + v*u.transpose());
  MATRIX3 Pi5uPF = 2*F*u*u.transpose();
  MATRIX3 Pi5vPF = 2*F*v*v.transpose();

  VECTOR9 g6 = flatten(Pi6PF);
  VECTOR9 g5u = flatten(Pi5uPF);
  VECTOR9 g5v = flatten(Pi5vPF);

  MATRIX3 uOut = 2*u*u.transpose();
  MATRIX3 vOut = 2*v*v.transpose();
  MATRIX3 uvOut = u*v.transpose() + v*u.transpose(); 
  MATRIX9 Hi6 = kronIdentity(uvOut);
  MATRIX9 H5u = kronIdentity(uOut);
  MATRIX9 H5v = kronIdentity(vOut);

  REAL i3 = invariant3(F);
  if(i3<_zeroGuard && i3 > -_zeroGuard) return MATRIX9::Zero(); //return 0 if NaNs are approaching.
  REAL sgn = i3 > 0 ? -1.0 : 1.0;
  if(i6*i6/(i5u*i5v) - 1 < _zeroGuard && i6*i6/(i5u*i5v) - 1 > -_zeroGuard) return MATRIX9::Zero(); //return 0 if NaNs approaching 

  REAL sqrtuv = sqrt(i5u*i5v);  
  REAL sqrt1Diff = sqrt(1 - i6*i6/(i5u*i5v));
  REAL diff = sgn*acos(i6/sqrtuv) - _theta0;
  
  //See Mathematica for these derivatives these are looooong
  REAL DPsiDi6 = -2*_mu*sgn*diff/(sqrtuv*sqrt1Diff);
  REAL DPsiDi5u = i5v*i6*_mu*sgn*diff/(sqrtuv*sqrtuv*sqrtuv*sqrt1Diff);
  REAL DPsiDi5v = i5u*i6*_mu*sgn*diff/(sqrtuv*sqrtuv*sqrtuv*sqrt1Diff);

  REAL D2PsiDi62 = 2*_mu*sgn*(i5u*i5v*sgn/(i5u*i5v-i6*i6) - i6*diff/(sqrtuv*sqrt1Diff*sqrt1Diff*sqrt1Diff))/(i5u*i5v);
  REAL D2PsiDi5u2 = i6*_mu*sgn*(sqrtuv*i6*sqrt1Diff*sgn+3*i5u*i5v*_theta0-2*i6*i6*_theta0 + (-3*i5u*i5v*sgn+2*i6*i6*sgn)*acos(i6/sqrtuv))/(2*i5u*i5u*sqrtuv*sqrt1Diff*(i5u*i5v-i6*i6));
  REAL D2PsiDi5v2 = i6*_mu*sgn*(sqrtuv*i6*sqrt1Diff*sgn+3*i5u*i5v*_theta0-2*i6*i6*_theta0 + (-3*i5u*i5v*sgn+2*i6*i6*sgn)*acos(i6/sqrtuv))/(2*i5v*i5v*sqrtuv*sqrt1Diff*(i5u*i5v-i6*i6));
  
  REAL D2PsiDi6i5u = -i5v*_mu*sgn*(sqrtuv*i6*sqrt1Diff*sgn + i5u*i5v*_theta0 - i5u*i5v*sgn*acos(i6/sqrtuv))/(sqrtuv*sqrtuv*sqrtuv*sqrt1Diff*(i5u*i5v-i6*i6));
  REAL D2PsiDi6i5v = -i5u*_mu*sgn*(sqrtuv*i6*sqrt1Diff*sgn + i5u*i5v*_theta0 - i5u*i5v*sgn*acos(i6/sqrtuv))/(sqrtuv*sqrtuv*sqrtuv*sqrt1Diff*(i5u*i5v-i6*i6));
  REAL D2PsiDi5ui5v = i6*_mu*sgn*(sqrtuv*i6*sqrt1Diff*sgn + i5u*i5v*_theta0 - i5u*i5v*sgn*acos(i6/sqrtuv))/(2*sqrtuv*sqrtuv*sqrtuv*sqrt1Diff*(i5u*i5v-i6*i6));

  return DPsiDi6*Hi6 + DPsiDi5u*H5u + DPsiDi5v*H5v
         + D2PsiDi62 * g6*g6.transpose() + D2PsiDi5u2 * g5u*g5u.transpose() + D2PsiDi5v2 * g5v*g5v.transpose()
         + D2PsiDi6i5u * (g6*g5u.transpose() + g5u * g6.transpose())
         + D2PsiDi6i5v * (g6*g5v.transpose() + g5v * g6.transpose())
         + D2PsiDi5ui5v * (g5u*g5v.transpose() + g5v*g5u.transpose());
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F (clamped)
// Analytical eigenfiltering done using slightly modified version
// of Haomiao's bending analysis
///////////////////////////////////////////////////////////////////////
MATRIX9 TET_STRAND_TWIST::clampedHessian(const MATRIX3& F) const
{
  REAL i3 = invariant3(F);
  if(i3<_zeroGuard && i3 > -_zeroGuard) return MATRIX9::Zero(); //return 0 if NaNs are approaching.
  REAL sgn = i3 > 0 ? -1.0 : 1.0;

  const VECTOR3 e0 = F.col(0);
  const VECTOR3 e1 = F.col(1);

  const REAL i6 = e0.dot(e1);
  const REAL i5u = e0.squaredNorm();
  const REAL i5v = e1.squaredNorm();
  if(i6*i6/(i5u*i5v) - 1 < _zeroGuard && i6*i6/(i5u*i5v) - 1 > -_zeroGuard) return MATRIX9::Zero(); //return 0 if NaNs approaching 
  
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const VECTOR3 crs = e0.cross(e1);
  REAL cnorm = crs.norm();
  const VECTOR3 n = crs / cnorm;
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  REAL dot0 = e0.dot(e0), dot1 = e1.dot(e1);
  REAL norm0 = e0.norm(), norm1 = e1.norm();
  //const REAL cosine = cos(theta * 0.5), sine = sin(theta * 0.5), tangent = tan(theta * 0.5);

  // partial psi partial theta (negated from introduction of strange sign ambiguities stuff I need to look into)
  const REAL g = 2*_mu * (sgn*theta - _theta0);
  const REAL h = 2*_mu * sgn;
  // cout << "g: "<<g<<"h: "<<h<<endl;
  // the two eigenvectors composed of n
  const REAL beta = (norm1/norm0 - norm0/norm1) * cos(theta);
  const REAL alpha = (- beta + sqrt(beta * beta + 4.0)) * 0.5;
  VECTOR9 q4;
  q4 << alpha * n,
        n,
        VECTOR3::Zero();
  const REAL lambda4 = g*(1.0/dot1/tan(theta) - alpha/norm0/norm1/sin(theta));
  VECTOR9 q5;
  q5 << n, 
       -alpha * n,
       VECTOR3::Zero();
  const REAL lambda5 = g*(1.0/dot1/tan(theta) + 1.0/alpha/norm0/norm1/sin(theta));
  // the rest 4
  const REAL n0 = h / dot0, n1 = h / dot1, w = g / h;
  const REAL w2 = w * w, nDiff = n0 - n1, nSum = n0 + n1;
  REAL gamma = nDiff/nSum; gamma = gamma * gamma;
  const REAL r = sqrt(4.0 * w2 * gamma + 1.0) * nSum;
  const REAL Rminus = sqrt(2.0 * (2.0 * w2 + 1.0 - r/nSum)) * nSum;
  const REAL Rplus = sqrt(2.0 * (2.0 * w2 + 1.0 + r/nSum)) * nSum;
  // cout << "r: "<<r<<"Rminus: "<<Rminus<<"nSum: "<<nSum<<endl;
  const REAL lambda0 = 0.25 * (-r - Rminus + nSum);
  const REAL lambda1 = 0.25 * (-r + Rminus + nSum);
  const REAL lambda2 = 0.25 * (r - Rplus + nSum);
  const REAL lambda3 = 0.25 * (r + Rplus + nSum);
  VECTOR9 Lam;
  Lam.setZero();
  Lam[0] = lambda0; Lam[1] = lambda1; Lam[2] = lambda2; Lam[3] = lambda3; Lam[4] = lambda4; Lam[5] = lambda5;
  Lam = sgn*Lam;
  for (int i = 0; i < 6; i++){
    Lam[i] = Lam[i] < 0? 0 : Lam[i];
  }
  VECTOR4 q0, q1, q2, q3;
  // q0.setZero(); q1.setZero(); q2.setZero(); q3.setZero(); 
  q0[0] = w/n1*(nDiff * nDiff - nSum * r - nDiff * Rminus);
  q0[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r+0.5*(nDiff-r)*Rminus)/n1;
  q0[2] = 4.0 * n1 * w;
  q0[3] = nSum - r - Rminus;
  q1[0] = w/n1*(nDiff * nDiff - nSum * r + nDiff * Rminus);
  q1[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r-0.5*(nDiff-r)*Rminus)/n1;
  q1[2] = 4.0 * n1 * w;
  q1[3] = nSum - r + Rminus;
  q2[0] = w/n1*(nDiff * nDiff + nSum * r - nDiff * Rplus);
  q2[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r+0.5*(nDiff+r)*Rplus)/n1;
  q2[2] = 4.0 * n1 * w;
  q2[3] = nSum + r - Rplus;
  q3[0] = w/n1*(nDiff * nDiff + nSum * r + nDiff * Rplus);
  q3[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r-0.5*(nDiff+r)*Rplus)/n1;
  q3[2] = 4.0 * n1 * w;
  q3[3] = nSum + r + Rplus;
  MATRIX9x4 A;
  A << e0, e0perp, VECTOR3::Zero(), VECTOR3::Zero(),
        VECTOR3::Zero(), VECTOR3::Zero(), e1, e1perp,
        VECTOR3::Zero(), VECTOR3::Zero(), VECTOR3::Zero(), VECTOR3::Zero();
  VECTOR9 qq0 = A * q0, qq1 = A * q1, qq2 = A * q2, qq3 = A * q3;
  qq0 = qq0.normalized(); qq1 = qq1.normalized(); qq2 = qq2.normalized(); 
  qq3 = qq3.normalized(); q4 = q4.normalized(); q5 = q5.normalized();
  MATRIX9 Q;
  Q.setZero();
  Q.col(0) = qq0; Q.col(1) = qq1; Q.col(2) = qq2; Q.col(3) = qq3; Q.col(4) = q4; Q.col(5) = q5;
  Q(6,6) = 1; Q(7,7) = 1; Q(8,8) = 1;
  MATRIX9 filtered = Q*Lam.asDiagonal()*Q.transpose();

  return filtered;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F (clamped)
// Analytical eigenfiltering done using slightly modified version
// of Haomiao's bending analysis
///////////////////////////////////////////////////////////////////////
void TET_STRAND_TWIST::clampedHessian(const MATRIX3& F, MATRIX9& plusH, MATRIX9& minusH) const
{
  REAL i3 = invariant3(F);
  if(i3<_zeroGuard && i3 > -_zeroGuard) 
  {
    plusH = minusH = MATRIX9::Zero();
    return; //return 0 if NaNs are approaching.
  } 
  REAL sgn = i3 > 0 ? -1.0 : 1.0;

  const VECTOR3 e0 = F.col(0);
  const VECTOR3 e1 = F.col(1);

  const REAL i6 = e0.dot(e1);
  const REAL i5u = e0.squaredNorm();
  const REAL i5v = e1.squaredNorm();
  if(i6*i6/(i5u*i5v) - 1 < _zeroGuard && i6*i6/(i5u*i5v) - 1 > -_zeroGuard) 
  {
    plusH = minusH = MATRIX9::Zero();
    return; //return 0 if NaNs are approaching.
  }

  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const VECTOR3 crs = e0.cross(e1);
  REAL cnorm = crs.norm();
  const VECTOR3 n = crs / cnorm;
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  REAL dot0 = e0.dot(e0), dot1 = e1.dot(e1);
  REAL norm0 = e0.norm(), norm1 = e1.norm();
  //const REAL cosine = cos(theta * 0.5), sine = sin(theta * 0.5), tangent = tan(theta * 0.5);

  // partial psi partial theta (negated from introduction of strange sign ambiguities stuff I need to look into)
  // TK: pushed the 2* _mu coefficient to later multiplies, because putting it here
  //       prevents you from setting _mu = 0. Computing w = g/h later then causes a divide-by-zero.
  const REAL g = (sgn*theta - _theta0);
  const REAL h = sgn;
  // cout << "g: "<<g<<"h: "<<h<<endl;
  // the two eigenvectors composed of n
  const REAL beta = (norm1/norm0 - norm0/norm1) * cos(theta);
  const REAL alpha = (- beta + sqrt(beta * beta + 4.0)) * 0.5;
  VECTOR9 q4;
  q4 << alpha * n,
        n,
        VECTOR3::Zero();
  const REAL lambda4 = 2 * _mu * g *(1.0/dot1/tan(theta) - alpha/norm0/norm1/sin(theta));
  VECTOR9 q5;
  q5 << n, 
       -alpha * n,
       VECTOR3::Zero();
  const REAL lambda5 = 2 * _mu * g *(1.0/dot1/tan(theta) + 1.0/alpha/norm0/norm1/sin(theta));
  // the rest 4
  const REAL hd0 = h / dot0;
  const REAL hd1 = h / dot1;
  const REAL n0 = 2 * _mu * hd0;
  const REAL n1 = 2 * _mu * hd1;
  const REAL w = g / h;
  const REAL w2 = w * w;
  const REAL nDiff = n0 - n1;
  const REAL nSum = n0 + n1;

  //REAL gamma = nDiff/nSum; 
  REAL gamma = (h / dot0 - h / dot1) / (h/dot0 + h / dot1);
  gamma = gamma * gamma;
  const REAL rad = sqrt(4.0 * w2 * gamma + 1.0);
  const REAL r = rad * nSum;
  const REAL Rminus = sqrt(2.0 * (2.0 * w2 + 1.0 - rad)) * nSum;
  const REAL Rplus = sqrt(2.0 * (2.0 * w2 + 1.0 + rad)) * nSum;
  // cout << "r: "<<r<<"Rminus: "<<Rminus<<"nSum: "<<nSum<<endl;
  const REAL lambda0 = 0.25 * (-r - Rminus + nSum);
  const REAL lambda1 = 0.25 * (-r + Rminus + nSum);
  const REAL lambda2 = 0.25 * (r - Rplus + nSum);
  const REAL lambda3 = 0.25 * (r + Rplus + nSum);
  VECTOR9 Lam;
  Lam.setZero();
  Lam[0] = lambda0; Lam[1] = lambda1; Lam[2] = lambda2; Lam[3] = lambda3; Lam[4] = lambda4; Lam[5] = lambda5;
  Lam = sgn*Lam;

  VECTOR9 LamPlus = Lam;
  VECTOR9 LamMinus = Lam;
  for (int i = 0; i < 6; i++){
    LamPlus[i] = LamPlus[i] < 0? 0 : LamPlus[i];
    LamMinus[i] = LamMinus[i] > 0? 0 : LamMinus[i];
  }
  VECTOR4 q0, q1, q2, q3;

  // TK: this one generates a NaN when mu = 0
#if 0
  q0[0] = w/n1*(nDiff * nDiff - nSum * r - nDiff * Rminus);
  q0[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r+0.5*(nDiff-r)*Rminus)/n1;
  q0[2] = 4.0 * n1 * w;
  q0[3] = nSum - r - Rminus;
  q1[0] = w/n1*(nDiff * nDiff - nSum * r + nDiff * Rminus);
  q1[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r-0.5*(nDiff-r)*Rminus)/n1;
  q1[2] = 4.0 * n1 * w;
  q1[3] = nSum - r + Rminus;
  q2[0] = w/n1*(nDiff * nDiff + nSum * r - nDiff * Rplus);
  q2[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r+0.5*(nDiff+r)*Rplus)/n1;
  q2[2] = 4.0 * n1 * w;
  q2[3] = nSum + r - Rplus;
  q3[0] = w/n1*(nDiff * nDiff + nSum * r + nDiff * Rplus);
  q3[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r-0.5*(nDiff+r)*Rplus)/n1;
  q3[2] = 4.0 * n1 * w;
  q3[3] = nSum + r + Rplus;
#else
  const REAL twoMuSq = 4.0 * _mu * _mu;

  q0[0] = w/hd1*(nDiff * nDiff - nSum * r - nDiff * Rminus);
  q0[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r+0.5*(nDiff-r)*Rminus)/hd1;
  q0[2] = (4.0 * hd1 * w) * twoMuSq;
  q0[3] = (nSum - r - Rminus) * (2.0 * _mu);

  q1[0] = w/hd1*(nDiff * nDiff - nSum * r + nDiff * Rminus);
  q1[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r-0.5*(nDiff-r)*Rminus)/hd1;
  q1[2] = (4.0 * hd1 * w) * twoMuSq;
  q1[3] = (nSum - r + Rminus) * (2.0 * _mu);

  q1[0] = w/hd1*(nDiff * nDiff - nSum * r + nDiff * Rminus);
  q1[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r-0.5*(nDiff-r)*Rminus)/hd1;
  q1[2] = (4.0 * hd1 * w) * twoMuSq;
  q1[3] = (nSum - r + Rminus) * (2.0 * _mu);
  
  q2[0] = w/hd1*(nDiff * nDiff + nSum * r - nDiff * Rplus);
  q2[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r+0.5*(nDiff+r)*Rplus)/hd1;
  q2[2] = (4.0 * hd1 * w) * twoMuSq;
  q2[3] = (nSum + r - Rplus) * (2.0 * _mu);

  q3[0] = w/hd1*(nDiff * nDiff + nSum * r + nDiff * Rplus);
  q3[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r-0.5*(nDiff+r)*Rplus)/hd1;
  q3[2] = (4.0 * hd1 * w) * twoMuSq;
  q3[3] = (nSum + r + Rplus) * (2.0 * _mu);

#endif

  MATRIX9x4 A;
  A << e0, e0perp, VECTOR3::Zero(), VECTOR3::Zero(),
        VECTOR3::Zero(), VECTOR3::Zero(), e1, e1perp,
        VECTOR3::Zero(), VECTOR3::Zero(), VECTOR3::Zero(), VECTOR3::Zero();
  VECTOR9 qq0 = A * q0, qq1 = A * q1, qq2 = A * q2, qq3 = A * q3;
  qq0 = qq0.normalized(); qq1 = qq1.normalized(); qq2 = qq2.normalized(); 
  qq3 = qq3.normalized(); q4 = q4.normalized(); q5 = q5.normalized();
  MATRIX9 Q;
  Q.setZero();
  Q.col(0) = qq0; Q.col(1) = qq1; Q.col(2) = qq2; Q.col(3) = qq3; Q.col(4) = q4; Q.col(5) = q5;
  Q(6,6) = 1; Q(7,7) = 1; Q(8,8) = 1;
  plusH  = Q*LamPlus.asDiagonal()*Q.transpose();
  minusH = Q*LamMinus.asDiagonal()*Q.transpose();
}

bool TET_STRAND_TWIST::energyNeedsSVD() const
{
  return true;
}

bool TET_STRAND_TWIST::PK1NeedsSVD() const
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX3 TET_STRAND_TWIST::computeFtwist(const vector<VECTOR3>& vertices)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  VECTOR3 e0 = vertices[0] - vertices[1];
  VECTOR3 e2 = vertices[3] - vertices[2];
  REAL alpha, beta;
  computeAlphaBeta(vertices, alpha, beta);

  e0 = e0 - alpha * e1;
  e2 = e2 - beta * e1;
  MATRIX3 Ftwist;
  Ftwist.col(0) = e0;
  Ftwist.col(1) = e2;
  Ftwist.col(2) = e1;

  return Ftwist;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 TET_STRAND_TWIST::computeTwistPFPx(const vector<VECTOR3>& vertices)
{
  REAL alpha, beta;
  computeAlphaBeta(vertices, alpha, beta);

  MATRIX9x12 pFpX;
  pFpX.setZero();
  pFpX(0,0) = 1;
  pFpX(0,3) = -1 + alpha;
  pFpX(0,6) = -alpha;
  pFpX(1,1) = 1;
  pFpX(1,4) = -1 + alpha;
  pFpX(1,7) = -alpha;
  pFpX(2,2) = 1;
  pFpX(2,5) = -1 + alpha;
  pFpX(2,8) = -alpha;
  pFpX(3,3) = beta;
  pFpX(3,6) = -1 - beta;
  pFpX(3,9) = 1;
  pFpX(4,4) = beta;
  pFpX(4,7) = -1 - beta;
  pFpX(4,10) = 1;
  pFpX(5,5) = beta;
  pFpX(5,8) = -1 - beta;
  pFpX(5,11) = 1;
  pFpX(6,3) = -1;
  pFpX(6,6) = 1;
  pFpX(7,4) = -1;
  pFpX(7,7) = 1;
  pFpX(8,5) = -1;
  pFpX(8,8) = 1;

  return pFpX;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TET_STRAND_TWIST::computeAlphaBeta(const vector<VECTOR3>& vertices, REAL& alpha, REAL& beta)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
    
  alpha = e0.dot(e1)/e1.squaredNorm();
  beta = e2.dot(e1)/e1.squaredNorm();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 TET_STRAND_TWIST::computeDelTwistPFPx(const vector<VECTOR3>& vertices)
{
  REAL alpha,beta;
  computeAlphaBeta(vertices, alpha, beta);
  VECTOR3 e0 = vertices[0]-vertices[1];
  VECTOR3 e1 = vertices[2]-vertices[1];
  VECTOR3 e2 = vertices[3]-vertices[2];

  // all the nontrivial del term blocks (with 1/e1^2 scalar factored out)
  MATRIX3 De0px0 = -e1*e1.transpose(); ///also De2px3
  MATRIX3 De0px1 = -e1*((2*alpha - 1)*e1 - e0).transpose();
  MATRIX3 De0px2 = -e1*(e0 - 2*alpha*e1).transpose();
  MATRIX3 De2px1 = -e1*(2*beta*e1 - e2).transpose();
  MATRIX3 De2px2 = e1*((2*beta + 1)*e1 - e2).transpose();

  MATRIX9x12 delpFpX;
  delpFpX.setZero();
  delpFpX.block<3,3>(0,0) = De0px0;
  delpFpX.block<3,3>(0,3) = De0px1;
  delpFpX.block<3,3>(0,6) = De0px2;
  delpFpX.block<3,3>(3,3) = De2px1;
  delpFpX.block<3,3>(3,6) = De2px2;
  delpFpX.block<3,3>(3,9) = De0px0;

  return 1.0/e1.squaredNorm()*delpFpX;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TET_STRAND_TWIST::computeFandGradients(const vector<VECTOR3>& vertices, 
                                            MATRIX3& Ftwist, 
                                            MATRIX9x12& pFpX, 
                                            MATRIX9x12& delpFpX)
{
  assert(vertices.size() == 4);
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
    
  const REAL alpha = e0.dot(e1)/e1.squaredNorm();
  const REAL beta  = e2.dot(e1)/e1.squaredNorm();
  
  // let's get F for twist
  const VECTOR3 e0perp = e0 - alpha * e1;
  const VECTOR3 e2perp = e2 - beta * e1;
  Ftwist.col(0) = e0perp;
  Ftwist.col(1) = e2perp;
  Ftwist.col(2) = e1;
  
  // let's get the projected gradient
  pFpX.setZero();
  pFpX(0,0) = 1.0;
  pFpX(0,3) = -1.0 + alpha;
  pFpX(0,6) = -alpha;
  pFpX(1,1) = 1.0;
  pFpX(1,4) = -1.0 + alpha;
  pFpX(1,7) = -alpha;
  pFpX(2,2) = 1.0;
  pFpX(2,5) = -1.0 + alpha;
  pFpX(2,8) = -alpha;
  pFpX(3,3) = beta;
  pFpX(3,6) = -1.0 - beta;
  pFpX(3,9) = 1.0;
  pFpX(4,4) = beta;
  pFpX(4,7) = -1.0 - beta;
  pFpX(4,10) = 1.0;
  pFpX(5,5) = beta;
  pFpX(5,8) = -1.0 - beta;
  pFpX(5,11) = 1.0;
  pFpX(6,3) = -1.0;
  pFpX(6,6) = 1.0;
  pFpX(7,4) = -1.0;
  pFpX(7,7) = 1.0;
  pFpX(8,5) = -1.0;
  pFpX(8,8) = 1.0;

  // let's get the parallel gradient
  const MATRIX3 De0px0 = -e1*e1.transpose(); ///also De2px3
  const MATRIX3 De0px1 = -e1*((2.0*alpha - 1.0)*e1 - e0).transpose();
  const MATRIX3 De0px2 = -e1*(e0 - 2.0*alpha*e1).transpose();
  const MATRIX3 De2px1 = -e1*(2.0*beta*e1 - e2).transpose();
  const MATRIX3 De2px2 = e1*((2.0*beta + 1.0)*e1 - e2).transpose();

  delpFpX.setZero();
  delpFpX.block<3,3>(0,0) = De0px0;
  delpFpX.block<3,3>(0,3) = De0px1;
  delpFpX.block<3,3>(0,6) = De0px2;
  delpFpX.block<3,3>(3,3) = De2px1;
  delpFpX.block<3,3>(3,6) = De2px2;
  delpFpX.block<3,3>(3,9) = De0px0;

  delpFpX *= 1.0/e1.squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// spatial Hessian, contains significant subtlety for twist
///////////////////////////////////////////////////////////////////////
MATRIX12 TET_STRAND_TWIST::forceGradient(const vector<VECTOR3>& vertices) const
{
  assert(vertices.size() == 4);

  //const MATRIX3 F = computeFtwist(vertices);
  //const MATRIX9x12 pFpX = computeTwistPFPx(vertices);
  //const MATRIX9x12 delpFpX = computeDelTwistPFPx(vertices);
  MATRIX3 F;
  MATRIX9x12 pFpX;
  MATRIX9x12 delpFpX;
  computeFandGradients(vertices, F, pFpX, delpFpX);

  const MATRIX9 H = hessian(F);
  return pFpX.transpose() * H * pFpX - delpFpX.transpose() * H * delpFpX;
}

///////////////////////////////////////////////////////////////////////
// spatial Hessian, contains significant subtlety for twist
///////////////////////////////////////////////////////////////////////
MATRIX12 TET_STRAND_TWIST::forceGradientRank4(const vector<VECTOR3>& vertices) const
{
  assert(vertices.size() == 4);

  const MATRIX3 F = computeFtwist(vertices);
  //const MATRIX9x12 pFpX = computeTwistPFPx(vertices);
  //const MATRIX9x12 delpFpX = computeDelTwistPFPx(vertices);
  const MATRIX9x12 pFpX = computeTwistPFPxFull(vertices);

  const MATRIX9 H = hessian(F);
  return pFpX.transpose() * H * pFpX + computeRank4Correction(vertices);
}

//////////////////////////////////////////////////////////////////////////////
// unsimplified original PFPx
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 TET_STRAND_TWIST::computeTwistPFPxFull(const vector<VECTOR3>& vertices)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
  const VECTOR3 e0perp = e0 - e0.dot(e1)/e1.squaredNorm()*e1;
  const VECTOR3 e2perp = e2 - e2.dot(e1)/e1.squaredNorm()*e1;
  const REAL norm1 = e1.norm();
  const VECTOR3 t1 = e1/norm1;
  const MATRIX3 z33 = MATRIX3::Zero(), I3 = MATRIX3::Identity();
  const MATRIX3 sigma = I3 - t1*t1.transpose();
  const MATRIX3 eta0 = -(t1.dot(e0)*sigma+t1*e0perp.transpose())/norm1;
  const MATRIX3 eta2 = -(t1.dot(e2)*sigma+t1*e2perp.transpose())/norm1;

  MATRIX9x12 pFpX;
  pFpX << sigma, -eta0-sigma, eta0, z33,
          z33, -eta2, eta2 - sigma, sigma,
          z33, -I3, I3, z33;

  return pFpX;
}

///////////////////////////////////////////////////////////////////////
// spatial Hessian, contains significant subtlety for twist
///////////////////////////////////////////////////////////////////////
MATRIX12 TET_STRAND_TWIST::clampedForceGradient(const vector<VECTOR3>& vertices) const
{
  assert(vertices.size() == 4);
  MATRIX3 F;
  MATRIX9x12 pFpX;
  MATRIX9x12 delpFpX;
  computeFandGradients(vertices, F, pFpX, delpFpX);

  MATRIX9 plusH, minusH;
  clampedHessian(F, plusH, minusH);
  //return pFpX.transpose() * plusH * pFpX - delpFpX.transpose() * minusH * delpFpX;

  const MATRIX12 finalH = pFpX.transpose() * plusH * pFpX - delpFpX.transpose() * minusH * delpFpX;

  if (finalH.hasNaN())
  {
    cout << " TET_STRAND_TWIST has NaN: " << endl;
    cout << " pFpX: " << endl << pFpX << endl;
    cout << " F: " << endl << F << endl;
    cout << " delpFpX: " << endl << delpFpX << endl;
    cout << " plusH: " << endl << plusH << endl;
    cout << " minusH: " << endl << minusH << endl;

    exit(0);
  }

  return finalH;
}

///////////////////////////////////////////////////////////////////////
// spatial Hessian, contains significant subtlety for twist
///////////////////////////////////////////////////////////////////////
MATRIX12 TET_STRAND_TWIST::clampedForceGradientRank4(const vector<VECTOR3>& vertices) const
{
  assert(vertices.size() == 4);
  const MATRIX3 F = computeFtwist(vertices);
  const MATRIX9x12 pFpX = computeTwistPFPxFull(vertices);
  //const MATRIX9 H = clampEigenvalues(hessian(F));
  //return pFpX.transpose() * H * pFpX + clampEigenvalues(computeRank4CorrectionClamped(vertices));
  const MATRIX9 H = clampedHessian(F);
  return pFpX.transpose() * H * pFpX + computeRank4CorrectionClamped(vertices);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 TET_STRAND_TWIST::computeRank4CorrectionClamped(const vector<VECTOR3>& vertices) const
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
  const VECTOR3 z0 = e0 - e0.dot(e1)/e1.squaredNorm()*e1;
  const VECTOR3 z1 = e2 - e2.dot(e1)/e1.squaredNorm()*e1;
  const REAL norm1 = e1.norm();
  const VECTOR3 t1 = e1/norm1;
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const VECTOR3 tb = z0.cross(z1).normalized();
  const VECTOR3 tau0 = z0 / normz0, tau1 = z1 / normz1;
  //const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized(); 
  //const MATRIX3  z33 = MATRIX3::Zero(); 
  const VECTOR3 z3 = VECTOR3::Zero();

  const REAL alpha = -4.0 * (tau0.dot(tau1perp)) / norm1 / norm1;
  const REAL beta = 2.0 * (tau0.dot(tau1perp))/norm1; 
  const REAL beta2 = beta*beta, alpha2 = alpha*alpha;
  const REAL a1 = t1.dot(e2) / norm1/normz1, b1 = 1/normz1;
  const REAL a0 = t1.dot(e0) / norm1/normz0, b0 = 1/normz0;
  VECTOR12 p0, p1, v0, v1, u0, u1;
  p0 << -b0*t1, (-a1-a0+b0)*t1, (a1+a0+b1)*t1, -b1*t1;
  p1 << b0*t1, (-a1+a0-b0)*t1, (a1-a0+b1)*t1, -b1*t1;
  v0 << z3, tau0, -tau0, z3; v1 << z3, tau1, -tau1, z3; 
  u0 = v0 + v1; u1 = v0 - v1;
  const REAL d00 = (p0.dot(p0)) / (u0.dot(u0)), d01 = (p0.dot(p1)) / (u1.dot(u1));
  const REAL d10 = (p0.dot(p1)) / (u0.dot(u0)), d11 = (p1.dot(p1)) / (u1.dot(u1));
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  SolveP4De(x, -alpha2 - beta2*(d11+d00), beta2*alpha*(d00-d11), beta2*beta2*(d11*d00 - d10*d01));
  MATRIX12 clamped; clamped.setZero();
  // this computes p Psi p theta. should be computed inside the material class.
  const REAL s = tb.dot(t1) > 0?-1:1;
  const REAL thetaGrad = 2.0 * s * _mu * (acos(tau0.dot(tau1)) - _theta0);
  for(int i = 0; i<4; i++){
    if(thetaGrad * x[i] > 0.0){
      const REAL f = beta/x[i], e = (1.0 /f/f - d00 + alpha/beta/f)/d10;
      const VECTOR12 eigVec = u0 * 1.0 + u1 * e + p0 * f + p1 * (e*f);
      clamped += thetaGrad * x[i] * eigVec * eigVec.transpose(); 
    }
  }

  return clamped;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 TET_STRAND_TWIST::computeRank4Correction(const vector<VECTOR3>& vertices) const
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
  const VECTOR3 z0 = e0 - e0.dot(e1)/e1.squaredNorm()*e1;
  const VECTOR3 z1 = e2 - e2.dot(e1)/e1.squaredNorm()*e1;
  const REAL norm1 = e1.norm();
  const VECTOR3 t1 = e1/norm1;
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const VECTOR3 tb = z0.cross(z1).normalized();
  const VECTOR3 tau0 = z0 / normz0, tau1 = z1 / normz1;
  const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized(); 
  const MATRIX3  z33 = MATRIX3::Zero(); const VECTOR3 z3 = VECTOR3::Zero();
  // blocks
  const MATRIX3 t1tau0perp = t1 * tau0perp.transpose(), t1tau1perp = t1 * tau1perp.transpose();
  MATRIX3x12 tau0pSigmapx, tau1pSigmapx, eta2px0, eta3px0, tau0peta2px, tau1peta3px;
  VECTOR12 eta2px1, eta3px1;
  tau0pSigmapx << z33, t1tau0perp, -t1tau0perp, z33;
  tau1pSigmapx << z33, t1tau1perp, -t1tau1perp, z33;
  tau0pSigmapx *= 1.0/norm1; tau1pSigmapx *= 1.0/norm1;
  const VECTOR3 reflect2 = e0 - 2.0 * t1 * (t1.dot(e0)), reflect3 = e2 - 2.0 * t1 * (t1.dot(e2)); 
  eta2px0 << z33, reflect2*tau0perp.transpose(), -reflect2*tau0perp.transpose(), z33;
  eta2px1 << e1, -reflect2 - e1, reflect2, z3;
  eta3px0 << z33, reflect3*tau1perp.transpose(), -reflect3*tau1perp.transpose(), z33;
  eta3px1 << z3, -reflect3, reflect3 - e1, e1;
  tau0peta2px = 1.0/norm1/norm1 *(eta2px0 - tau0perp*eta2px1.transpose());
  tau1peta3px = 1.0/norm1/norm1 *(eta3px0 - tau1perp*eta3px1.transpose());
  // this computes p Psi p theta. should be computed inside the material class.
  const REAL s = tb.dot(t1) > 0?-1:1;
  const REAL thetaGrad = 2.0 * s * _mu * (acos(tau0.dot(tau1)) - _theta0);
  MATRIX12 extra;
  extra << 1.0/normz0 * tau0pSigmapx,
           1.0/normz0*(-tau0peta2px-tau0pSigmapx) + 1.0/normz1 * tau1peta3px,
           1.0/normz0*tau0peta2px - 1.0/normz1*(tau1peta3px - tau1pSigmapx),
          -1.0/normz1*tau1pSigmapx;
  return thetaGrad * extra;
}

} // VOLUME
} // HOBAK
