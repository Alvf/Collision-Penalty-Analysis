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
// May 26, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "QUADRATIC_F_BENDING.h"
#include "util/MATRIX_UTIL.h"
#include "TIMER.h"
#include "ext/solvePoly/poly34.h"

#include <iostream>
using namespace std;

#define USING_UNFILTERED 0
#define USING_GAUSS_NEWTON_BENDING 0
#define USING_BRUTE_FORCE_CLAMP 0
#define USING_RANK4_CORRECTION 1

namespace HOBAK {
namespace SHELL {

QUADRATIC_F_BENDING::QUADRATIC_F_BENDING(const REAL& mu) :
  BENDING(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string QUADRATIC_F_BENDING::name() const
{ 
  return std::string("QUADRATIC_F_BENDING"); 
}

REAL clampValue(const REAL v, const REAL up, const REAL down)
{
  if(v > up)
    return up;
  else if(v < down)
    return down;
  else
    return v;
}

REAL QUADRATIC_F_BENDING::psi(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  // compute theta
  const VECTOR3 y1 = flap[2] - flap[1];
  const VECTOR3 y2 = flap[0] - flap[1];
  const VECTOR3 y3 = flap[3] - flap[1];
  const VECTOR3 t1 = y1.normalized();
  const VECTOR3 z0 = y2 - y2.dot(t1) * t1;
  const VECTOR3 z1 = y3 - y3.dot(t1) * t1;
  const REAL s = t1.dot(z0.cross(z1)) > 0?1.0:-1.0;
  const REAL theta = s * (M_PI -  acos(clampValue(z0.dot(z1) / z0.norm() / z1.norm(), 1.0, -1.0)));
  const REAL diff = theta - restTheta;

  return 0.5 * _mu * diff * diff;
}

VECTOR12 QUADRATIC_F_BENDING::gradient(const std::vector<VECTOR3>& flap, const REAL& restTheta) const
{
  // compute theta
  const VECTOR3 y1 = flap[2] - flap[1];
  const VECTOR3 y2 = flap[0] - flap[1];
  const VECTOR3 y3 = flap[3] - flap[1];
  const REAL norm1 = y1.norm();
  const VECTOR3 t1 = y1/norm1;
  const VECTOR3 z0 = y2 - y2.dot(t1) * t1;
  const VECTOR3 z1 = y3 - y3.dot(t1) * t1;
  VECTOR3 tb = z0.cross(z1); tb = tb.normalized();
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const VECTOR3 tau0 = z0/normz0, tau1 = z1/normz1;
  const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized();
  const REAL s = t1.dot(tb) > 0?1.0:-1.0;
  const REAL theta = s * (M_PI -  acos(clampValue(z0.dot(z1) / normz0 / normz1, 1.0, -1.0)));
  const REAL diff = theta - restTheta, g = _mu * diff;
  // cout<<"theta: "<<theta<<" rest theta: "<<restTheta<<endl;

  // pFpx simplified
  // const MATRIX3 sigma = I3;
  // const MATRIX3 eta2 = -(t1.dot(y2)*sigma)/norm1;
  // const MATRIX3 eta3 = -(t1.dot(y3)*sigma)/norm1;
  // MATRIX9x12 pFpX;
  // pFpX << sigma, -eta2-sigma, eta2, z33,
  //         z33, -eta3-sigma, eta3, sigma,
  // pPsipF

  // coeffs
  const VECTOR3 g0 = -g * s/ normz0 * tau0perp, g1 = g * s/normz1 * tau1perp;
  //   if(g0.hasNaN()){
  //   cout<<"gradient nan.\n";
  //   cout<<"theta: "<<theta<<endl;
  //   exit(0);
  // }
  const REAL eta2 = -t1.dot(y2)/norm1, eta3 = -t1.dot(y3)/norm1;
  VECTOR12 G;
  G << g0, (-eta2 - 1) * g0 + (-eta3 - 1) * g1, eta2 * g0 + eta3 * g1, g1;
  return G;
}

MATRIX3 getCrossMatrix(const VECTOR3 v)
{
  MATRIX3 M;
  M << 0, v[2], -v[1],
       -v[2], 0, v[0],
       v[1], -v[0], 0;
  return M;
}

MATRIX12 QUADRATIC_F_BENDING::hessian(const std::vector<VECTOR3>& flap, const REAL& restTheta) const
{
  // compute theta
  const VECTOR3 y1 = flap[2] - flap[1];
  const VECTOR3 y2 = flap[0] - flap[1];
  const VECTOR3 y3 = flap[3] - flap[1];
  const REAL norm1 = y1.norm();
  const VECTOR3 t1 = y1/norm1;
  const VECTOR3 z0 = y2 - y2.dot(t1) * t1;
  const VECTOR3 z1 = y3 - y3.dot(t1) * t1;
  VECTOR3 tb = z0.cross(z1); tb = tb.normalized();
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const REAL dot0 = normz0 * normz0, dot1 = normz1 * normz1;
  const VECTOR3 tau0 = z0/normz0, tau1 = z1/normz1;
  const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized();
  const REAL crs = (tau0.cross(tau1)).norm();
  const REAL s = t1.dot(tb) > 0?1.0:-1.0;
  const REAL theta = s * (M_PI -  acos(clampValue(z0.dot(z1) / z0.norm() / z1.norm(), 1.0, -1.0)));
  const REAL diff = theta - restTheta, g = _mu * diff, h = _mu;

  // using the rank-2 correction
  // pFpx simplified
  const MATRIX3 z33 = MATRIX3::Zero();
  // const MATRIX3 sigmaSimp = MATRIX3::Identity();
  // const MATRIX3 eta2Simp = -(t1.dot(y2)*sigmaSimp)/norm1;
  // const MATRIX3 eta3Simp = -(t1.dot(y3)*sigmaSimp)/norm1;
  // MATRIX6x12 pFpXSimp;
  // pFpXSimp << sigmaSimp, -eta2Simp-sigmaSimp, eta2Simp, z33,
  //         z33, -eta3Simp-sigmaSimp, eta3Simp, sigmaSimp;
  // // pFpx original - pFpx simplified
  // const MATRIX3 sigmaDiff = -t1 * t1.transpose();
  // const MATRIX3 eta2Diff = -(t1.dot(y2)*sigmaDiff + t1 * z0.transpose())/norm1;
  // const MATRIX3 eta3Diff = -(t1.dot(y3)*sigmaDiff + t1 * z1.transpose())/norm1;
  // MATRIX6x12 pFpXDiff;
  // pFpXDiff << sigmaDiff, -eta2Diff-sigmaDiff, eta2Diff, z33,
  //         z33, -eta3Diff-sigmaDiff, eta3Diff, sigmaDiff;
  const MATRIX3 sigma = MATRIX3::Identity()-t1 * t1.transpose();
  const MATRIX3 eta2 = -(t1.dot(y2)*sigma + t1 * z0.transpose())/norm1;
  const MATRIX3 eta3 = -(t1.dot(y3)*sigma + t1 * z1.transpose())/norm1;
  MATRIX6x12 pFpX;
  pFpX << sigma, -eta2-sigma, eta2, z33,
          z33, -eta3-sigma, eta3, sigma;

  // pPsipF
  // theta gradient
  VECTOR6 G;
  const VECTOR3 g0 = tau0perp / normz0, g1 = -tau1perp / normz1;
  G << g0, g1;
  // theta hessian
  const MATRIX3 Tb = getCrossMatrix(tb), T0 = getCrossMatrix(tau0), T1 = getCrossMatrix(tau1);
  MATRIX6 H;
  if(crs == 0){
    H.setZero();
  }
  else{
    MATRIX3 B00, B01, B11, B10;
    B00 = 1.0/normz0/normz0 *(Tb - 2.0 * tau0perp * tau0.transpose()
                  -(T0 * T1 + tau0perp*tau1perp.transpose())/crs);
    B01 = 1.0/crs/normz0/normz1 * (T0 * T0 + tau0perp*tau0perp.transpose());
    B11 = 1.0/normz1/normz1 *(-Tb + 2.0 * tau1perp * tau1.transpose()
                  -(T1 * T0 + tau1perp*tau0perp.transpose())/crs);
    B10 = 1.0/crs/normz0/normz1 * (T1 * T1 + tau1perp*tau1perp.transpose());
    H << B00, B01,
         B10, B11;
  }
  const MATRIX6 pPsipF = h * G * G.transpose() - s * g * H;

  // blocks
  const VECTOR3 z3 = VECTOR3::Zero();
  const MATRIX3 t1tau0perp = t1 * tau0perp.transpose(), t1tau1perp = t1 * tau1perp.transpose();
  MATRIX3x12 tau0pSigmapx, tau1pSigmapx, eta2px0, eta3px0, tau0peta2px, tau1peta3px;
  VECTOR12 eta2px1, eta3px1;
  tau0pSigmapx << z33, t1tau0perp, -t1tau0perp, z33;
  tau1pSigmapx << z33, t1tau1perp, -t1tau1perp, z33;
  tau0pSigmapx *= 1.0/norm1; tau1pSigmapx *= 1.0/norm1;
  const VECTOR3 reflect2 = y2 - 2.0 * t1 * (t1.dot(y2)), reflect3 = y3 - 2.0 * t1 * (t1.dot(y3)); 
  eta2px0 << z33, reflect2*tau0perp.transpose(), -reflect2*tau0perp.transpose(), z33;
  eta2px1 << y1, -reflect2 - y1, reflect2, z3;
  eta3px0 << z33, reflect3*tau1perp.transpose(), -reflect3*tau1perp.transpose(), z33;
  eta3px1 << z3, -reflect3 - y1, reflect3, y1;
  tau0peta2px = 1.0/norm1/norm1 *(eta2px0 - tau0perp*eta2px1.transpose());
  tau1peta3px = 1.0/norm1/norm1 *(eta3px0 - tau1perp*eta3px1.transpose());
  // this computes p Psi p theta. should be computed inside the material class.
  // const REAL s = tb.dot(t1) > 0?-1:1;
  // const REAL thetaGrad = 2.0 * s * material->mu() * (acos(tau0.dot(tau1)) - material->theta0());
  MATRIX12 extra;
  extra << 1.0/normz0 * tau0pSigmapx,
           1.0/normz0*(-tau0peta2px-tau0pSigmapx) + 1.0/normz1 * (tau1peta3px + tau1pSigmapx),
           1.0/normz0*tau0peta2px - 1.0/normz1*tau1peta3px,
          -1.0/normz1*tau1pSigmapx;
  // return thetaGrad * extra;
  MATRIX12 rank4 = pFpX.transpose() * pPsipF * pFpX - s * g * extra;

  // MATRIX12 rank2 = pFpXSimp.transpose() * pPsipF * pFpXSimp - pFpXDiff.transpose() * pPsipF * pFpXDiff;
  // cout<<"extra: \n"<<(- s * g *extra)<<endl;
  // cout<<"rank2 extra: \n"<< (rank2 - pFpX.transpose() * pPsipF * pFpX)<<endl;
  return rank4;
}

MATRIX12 QUADRATIC_F_BENDING::clampedHessian(const std::vector<VECTOR3>& flap, const REAL& restTheta) const
{
  TIMER functionTimer(__FUNCTION__);
  #if USING_BRUTE_FORCE_CLAMP
    return clampEigenvalues(hessian(flap, restTheta));
  #endif
  #if USING_UNFILTERED
    return hessian(flap, restTheta);
  #endif
  // compute theta
  const VECTOR3 y1 = flap[2] - flap[1];
  const VECTOR3 y2 = flap[0] - flap[1];
  const VECTOR3 y3 = flap[3] - flap[1];
  const REAL norm1 = y1.norm();
  const VECTOR3 t1 = y1/norm1;
  const VECTOR3 z0 = y2 - y2.dot(t1) * t1;
  const VECTOR3 z1 = y3 - y3.dot(t1) * t1;
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const REAL dot0 = normz0 * normz0, dot1 = normz1 * normz1;
  VECTOR3 tb = z0.cross(z1); tb = tb.normalized();
  const VECTOR3 tau0 = z0/normz0, tau1 = z1/normz1;
  const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized();
  const VECTOR3 z0perp = tau0perp * normz0, z1perp = tau1perp * normz1;
  const REAL crs = (tau0.cross(tau1)).norm();
  const REAL s = t1.dot(tb) > 0?1.0:-1.0;
  const REAL theta = s * (M_PI -  acos(clampValue(z0.dot(z1) / z0.norm() / z1.norm(), 1.0, -1.0)));
  const REAL diff = theta - restTheta, g = _mu * diff, h = _mu;

  const MATRIX3 z33 = MATRIX3::Zero();
  const MATRIX3 sigma = MATRIX3::Identity()-t1 * t1.transpose();
  const MATRIX3 eta2 = -(t1.dot(y2)*sigma + t1 * z0.transpose())/norm1;
  const MATRIX3 eta3 = -(t1.dot(y3)*sigma + t1 * z1.transpose())/norm1;
  MATRIX6x12 pFpX;
  pFpX << sigma, -eta2-sigma, eta2, z33,
          z33, -eta3-sigma, eta3, sigma;
  // const MATRIX3 sigmaSimp = MATRIX3::Identity();
  // const MATRIX3 eta2Simp = -(t1.dot(y2)*sigmaSimp)/norm1;
  // const MATRIX3 eta3Simp = -(t1.dot(y3)*sigmaSimp)/norm1;
  // MATRIX6x12 pFpXSimp;
  // pFpXSimp << sigmaSimp, -eta2Simp-sigmaSimp, eta2Simp, z33,
  //         z33, -eta3Simp-sigmaSimp, eta3Simp, sigmaSimp;
  #if USING_GAUSS_NEWTON_BENDING
    // pFpx simplified
    // pFpx original - pFpx simplified
    // const MATRIX3 sigmaDiff = -t1 * t1.transpose();
    // const MATRIX3 eta2Diff = -(t1.dot(y2)*sigmaDiff + t1 * z0.transpose())/norm1;
    // const MATRIX3 eta3Diff = -(t1.dot(y3)*sigmaDiff + t1 * z1.transpose())/norm1;
    // MATRIX6x12 pFpXDiff;
    // pFpXDiff << sigmaDiff, -eta2Diff-sigmaDiff, eta2Diff, z33,
    //         z33, -eta3Diff-sigmaDiff, eta3Diff, sigmaDiff;

    // pPsipF
    // theta gradient
    VECTOR6 G;
    const VECTOR3 g0 = tau0perp / normz0, g1 = -tau1perp / normz1;
    G << g0, g1;
    const MATRIX6 pPsipF = h * G * G.transpose();

    // return pFpXSimp.transpose() * pPsipF * pFpXSimp - pFpXDiff.transpose() * pPsipF * pFpXDiff;
    return pFpX.transpose() * pPsipF * pFpX;
  #endif

  // eigensystem of p2PsipF2
  // the two eigenvectors composed of tb
  VECTOR6 lambda = VECTOR6::Zero(), lambdaPos = VECTOR6::Zero(), lambdaNeg = VECTOR6::Zero();
  const REAL beta = -(normz1/normz0 - normz0/normz1) * cos(theta);
  const REAL alpha = (- beta + sqrt(beta * beta + 4.0)) * 0.5;
  VECTOR6 q4;
  q4 << alpha * tb, tb;
  VECTOR6 q5;
  q5 << tb, -alpha * tb;
  if(abs(theta) < 1e-16){
    lambda(4) = 0.0;
    lambda(5) = 0.0;
  }
  else{
    lambda(4) = -g/sin(theta)*(-cos(theta)/dot1 - alpha/normz0/normz1);
    lambda(5) = -g/sin(theta)*(-cos(theta)/dot1 + 1.0/alpha/normz0/normz1);
  }

  // the rest 4
  const REAL n0 = h / dot0, n1 = h / dot1, w = -s*g / h;
  const REAL w2 = w * w, nDiff = n0 - n1, nSum = n0 + n1;
  REAL gamma = nDiff/nSum; gamma = gamma * gamma;
  const REAL r = sqrt(4.0 * w2 * gamma + 1.0) * nSum;
  const REAL Rminus = sqrt(2.0 * (2.0 * w2 + 1.0 - r/nSum)) * nSum;
  const REAL Rplus = sqrt(2.0 * (2.0 * w2 + 1.0 + r/nSum)) * nSum;
// cout << "r: "<<r<<"Rminus: "<<Rminus<<"nSum: "<<nSum<<endl;
  lambda(0) = 0.25 * (-r - Rminus + nSum);
  lambda(1) = 0.25 * (-r + Rminus + nSum);
  lambda(2) = 0.25 * (r - Rplus + nSum);
  lambda(3) = 0.25 * (r + Rplus + nSum);

  VECTOR4 q0, q1, q2, q3;
  const VECTOR3 z3 = VECTOR3::Zero();
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
  MATRIX6x4 A;
  A << z0, z0perp, z3, z3,
        z3, z3, z1, z1perp;
  VECTOR6 qq0 = A * q0, qq1 = A * q1, qq2 = A * q2, qq3 = A * q3;
  qq0 = qq0.normalized(); qq1 = qq1.normalized(); qq2 = qq2.normalized(); 
  qq3 = qq3.normalized(); q4 = q4.normalized(); q5 = q5.normalized();
  MATRIX6 Q; Q << qq0, qq1, qq2, qq3, q4, q5;
  for(int i = 0; i < 6; i++) {
    if(lambda(i) > 0)
      lambdaPos(i) = lambda(i);
    else
      lambdaNeg(i) = lambda(i);
  } 

  const MATRIX6 filterPos = Q * lambdaPos.asDiagonal() * Q.transpose();
  // const MATRIX6 unfiltered = Q * lambda.asDiagonal() * Q.transpose();
  // cout<<"z0: "<<z0<<endl; cout<<"z1: "<<z1<<endl;
  // cout<<"lambda: "<<lambda<<endl;
  // cout<<"unfiltered: "<<unfiltered<<endl;
  // cout<<"Q: "<<Q<<endl;
  // using the rank 2 decomposition
  // pFpx simplified
  // const MATRIX3 z33 = MATRIX3::Zero();
  // const MATRIX3 sigmaSimp = MATRIX3::Identity();
  // const MATRIX3 eta2Simp = -(t1.dot(y2)*sigmaSimp)/norm1;
  // const MATRIX3 eta3Simp = -(t1.dot(y3)*sigmaSimp)/norm1;
  // MATRIX6x12 pFpXSimp;
  // pFpXSimp << sigmaSimp, -eta2Simp-sigmaSimp, eta2Simp, z33,
  //         z33, -eta3Simp-sigmaSimp, eta3Simp, sigmaSimp;
  // pFpx original - pFpx simplified
  // return pFpXSimp.transpose() * filterPos * pFpXSimp - pFpXDiff.transpose() * filterNeg * pFpXDiff;

  #if USING_RANK4_CORRECTION
  // TIMER rank4Timer("rank4 assembly");
    const REAL aa = -4.0 * (tau0.dot(tau1perp)) / norm1 / norm1;
    const REAL bb = 2.0 * (tau0.dot(tau1perp))/norm1; 
    const REAL bb2 = bb*bb, aa2 = aa*aa;
    const REAL a1 = t1.dot(y3) / norm1/normz1, b1 = 1/normz1;
    const REAL a0 = t1.dot(y2) / norm1/normz0, b0 = 1/normz0;
    VECTOR12 p0, p1, v0, v1, u0, u1;
    p0 << -b0*t1, (-a1-a0+b1+b0)*t1, (a1+a0)*t1, -b1*t1;
    p1 << b0*t1, (-a1+a0+b1-b0)*t1, (a1-a0)*t1, -b1*t1;
    v0 << z3, tau0, -tau0, z3; v1 << z3, tau1, -tau1, z3; 
    u0 = v0 + v1; u1 = v0 - v1;
    REAL d00, d01, d10, d11;
    const REAL dotu0 = u0.dot(u0), dotu1 = u1.dot(u1);
    if(dotu0 < 1e-22){
      d00 = 0.0; d10 = 0.0;
    }
    else{
      d00 = (p0.dot(p0)) / dotu0; d10 = (p0.dot(p1)) / dotu0;
    }
    if(dotu1 < 1e-22){
      d01 = 0.0; d11 = 0.0;
    }
    else{
      d01 = (p0.dot(p1)) / dotu1; d11 = (p1.dot(p1)) / dotu1;
    }
    // const REAL d00 = (p0.dot(p0)) / (u0.dot(u0)), d01 = (p0.dot(p1)) / (u1.dot(u1));
    // const REAL d10 = (p0.dot(p1)) / (u0.dot(u0)), d11 = (p1.dot(p1)) / (u1.dot(u1));
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    SolveP4De(x, -aa2 - bb2*(d11+d00), bb2*aa*(d00-d11), bb2*bb2*(d11*d00 - d10*d01));
    MATRIX12 clamped; clamped.setZero();
    // this computes p Psi p theta. should be computed inside the material class.
    const REAL thetaGrad = -s*g;
    REAL e, f;
    for(int i = 0; i<4; i++){
      if(x[i] * thetaGrad > 0.0){
        f = abs(x[i])<1e-16?0.0:bb/x[i];
        if(dotu0 < 1e-22 || abs(f)<1e-22)
          e = 0.0;
        else
          e = (1.0 /f/f - d00 + aa/bb/f)/d10;
        // const REAL f = bb/x[i], e = (1.0 /f/f - d00 + aa/bb/f)/d10;
        VECTOR12 eigVec = u0 * 1.0 + u1 * e + p0 * f + p1 * (e*f);
        eigVec= eigVec.normalized();
        clamped += thetaGrad * x[i] * eigVec * eigVec.transpose(); 
      }
    }
    // rank4Timer.stop();
    // MATRIX6x12 pFpX = pFpXSimp + pFpXDiff;
    // if(clamped.hasNaN()){
    //   cout<<"e: "<<e<<endl;
    //   cout<<"f: "<<f<<endl;
    //   cout<<"dotu0: "<<dotu0<<endl;
    //   cout<<"dotu1: "<<e<<endl;
    // }
    // MATRIX12 rank4 =  clamped + pFpX.transpose() * unfiltered * pFpX;
    return clamped + pFpX.transpose() * filterPos * pFpX;
  #endif
  // MATRIX12 rank2 = pFpXSimp.transpose() * unfiltered * pFpXSimp - pFpXDiff.transpose() * unfiltered * pFpXDiff;
  // MATRIX12 extraTerm = rank2 - pFpX.transpose() * unfiltered * pFpX;
  // Eigen::SelfAdjointEigenSolver<MATRIX12> eigensolver(extraTerm);
  // const MATRIX9 Q = eigensolver.eigenvectors();
  // VECTOR12 values = eigensolver.eigenvalues();
  // cout<<"extraterm eigenvalues: \n"<<values<<endl;
  // cout<<"clamped eigenvalues: \n"<<thetaGrad*x[0]<<" "<<thetaGrad*x[1]<<" "<<thetaGrad*x[2]<<" "<<thetaGrad*x[3]<<endl;
  // cout<<"extraterm: \n"<<extraTerm<<endl;
  // cout<<"clamped: \n"<<clamped<<endl;
  // exit(0);
  // return rank4;
  const MATRIX3 sigmaDiff = -t1 * t1.transpose();
  const MATRIX3 eta2Diff = -(t1.dot(y2)*sigmaDiff + t1 * z0.transpose())/norm1;
  const MATRIX3 eta3Diff = -(t1.dot(y3)*sigmaDiff + t1 * z1.transpose())/norm1;
  MATRIX6x12 pFpXDiff;
  pFpXDiff << sigmaDiff, -eta2Diff-sigmaDiff, eta2Diff, z33,
          z33, -eta3Diff-sigmaDiff, eta3Diff, sigmaDiff;
  const MATRIX6x12 pFpXSimp = pFpX - pFpXDiff;
  const MATRIX6 filterNeg = Q * lambdaNeg.asDiagonal() * Q.transpose();
  return pFpXSimp.transpose() * filterPos * pFpXSimp - pFpXDiff.transpose() * filterNeg * pFpXDiff;
}

REAL QUADRATIC_F_BENDING::restAngle(const std::vector<VECTOR3>& restFlap) const
{
  const VECTOR3 y1 = restFlap[2] - restFlap[1];
  const VECTOR3 y2 = restFlap[0] - restFlap[1];
  const VECTOR3 y3 = restFlap[3] - restFlap[1];
  const REAL norm1 = y1.norm();
  const VECTOR3 t1 = y1/norm1;
  const VECTOR3 z0 = y2 - y2.dot(t1) * t1;
  const VECTOR3 z1 = y3 - y3.dot(t1) * t1;
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  VECTOR3 tb = z0.cross(z1); tb = tb.normalized();
  const REAL s = t1.dot(tb) > 0?1.0:-1.0;
  const REAL theta = s * (M_PI -  acos(clampValue(z0.dot(z1) / z0.norm() / z1.norm(), 1.0, -1.0)));
  return theta;
}


} //SHELL
} //HOBAK