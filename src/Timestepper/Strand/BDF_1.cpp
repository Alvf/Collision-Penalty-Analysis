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
#include "BDF_1.h"
#include "TIMER.h"

using namespace std;

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
BDF_1::BDF_1(STRAND_MESH& strandMesh, 
             STRAND::STRETCHING& stretching,
             ENERGY_12D& edgeEdgeEnergy) :
  TIMESTEPPER(strandMesh, stretching, edgeEdgeEnergy)
{
  _accelerationOld.resize(_DOFs);
  _acceleration.resize(_DOFs);
  _accelerationOld.setZero();
  _acceleration.setZero();

  _maxNewtonIterations = 10;

  // these are the BDF-1 settings of the Newmark constants
  // this is here mostly for illustrative purposes. Many of the alphas
  // never get called, since we can assume they are zero.
	_alpha[0] = 1.0 / (_dt * _dt);
	_alpha[1] = 1.0 / _dt;
	_alpha[2] = 0;
	_alpha[3] = 1.0 / _dt;
	_alpha[4] = 0;
	_alpha[5] = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void BDF_1::setDt(const REAL dt)
{
  _dt = dt;
	_alpha[0] = 1.0 / (_dt * _dt);
	_alpha[1] = 1.0 / _dt;
	_alpha[2] = 0;
	_alpha[3] = 1.0 / _dt;
	_alpha[4] = 0;
	_alpha[5] = 0;
  cout << " New dt: " << _dt << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BDF_1::solveDynamics(const bool verbose)
{
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BDF-1 RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // store the internal forces for later
  VECTOR R;
  VECTOR z;

  // do Newton-Raphson
  REAL eps = _residualTolerance;
  REAL maxR = eps * 10;
  unsigned int step = 0;
  
  while (step < _maxNewtonIterations && maxR > eps)
  {
    // build new constraints and see if we should break any
    findNewSurfaceConstraints(verbose);
    buildConstraintMatrixFaster();

    // do collision detection, including spatial data structure updates 
    if (_collisionsEnabled)
      _strandMesh.computeEdgeEdgeCollisions(verbose);

    // "z is a vector of the desired values for the constrained variables",
    // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
    // _constraintTargets did not project off the null directions
    updateConstraintTargets();
    z =_IminusS * _constraintTargets;

    // get the internal forces
    TIMER forceTimer("Internal forces");
    const VECTOR stretchingForce = _strandMesh.computeStretchingForces(_stretchingEnergy);
    const VECTOR bendingForce = _strandMesh.computeBendingForces();
    const VECTOR twistingForce = _strandMesh.computeTwistingForces();

    if (verbose)
    {
      cout << " Bending force :   " << bendingForce.norm() << endl;
      cout << " Stretching force: " << stretchingForce.norm() << endl;
      cout << " Twisting force:   " << twistingForce.norm() << endl;
      cout << " Forming linear system ... " << flush;
    }

    // making non-const so we can add collision forces
    VECTOR R = bendingForce + stretchingForce + twistingForce;
    forceTimer.stop();

    // get the reduced stiffness matrix
    SPARSE_MATRIX K;
    if (_hessianClampingEnabled)
    {
      if (verbose)
      {
        cout << " USING CLAMPED HESSIAN" << endl;
        cout << " Building material Hessian ..." << flush;
      }
      K = _strandMesh.computeStretchingClampedHessian(_stretchingEnergy);
      K += _strandMesh.computeBendingClampedHessian();
      K += _strandMesh.computeTwistingClampedHessian();

      if (verbose)
        cout << " done." << endl;
    }
    else if (_gaussNewtonEnabled)
    {
      if (verbose)
      {
        cout << " USING GAUSS-NEWTON HESSIAN" << endl;
        cout << " Building material Hessian ..." << flush;
      }
      //K = _strandMesh.computeStretchingClampedHessian(_stretchingEnergy);
      K = _strandMesh.computeStretchingGaussNewtonHessian(_stretchingEnergy);
      //K += _strandMesh.computeBendingClampedHessian();
      K += _strandMesh.computeBendingGaussNewtonHessian();
      //K += _strandMesh.computeTwistingClampedHessian();
      K += _strandMesh.computeTwistingGaussNewtonHessian();

      if (verbose)
        cout << " done." << endl;
    }
    else if (_gaussNewtonFilteredEnabled)
    {
      if (verbose)
      {
        cout << " USING GAUSS-NEWTON HESSIAN, EXCEPT STRETCHING IS FILTERED" << endl;
        cout << " Building material Hessian ..." << flush;
      }
      K = _strandMesh.computeStretchingClampedHessian(_stretchingEnergy);
      //K += _strandMesh.computeBendingClampedHessian();
      K += _strandMesh.computeBendingGaussNewtonHessian();
      //K += _strandMesh.computeTwistingClampedHessian();
      K += _strandMesh.computeTwistingGaussNewtonHessian();

      if (verbose)
        cout << " done." << endl;
    }
    else
    {
      //SPARSE_MATRIX K = _strandMesh.computeStretchingHessian(_stretchingEnergy);
      //K += _strandMesh.computeBendingHessian();
      //K += _strandMesh.computeTwistingHessian();
      SPARSE_MATRIX Ks = _strandMesh.computeStretchingHessian(_stretchingEnergy);
      SPARSE_MATRIX Kb = _strandMesh.computeBendingHessian();
      SPARSE_MATRIX Kt = _strandMesh.computeTwistingHessian();
      if (verbose)
      {
        cout << " Bending hessian:    " << Kb.norm() << endl;
        cout << " Stretching hessian: " << Ks.norm() << endl;
        cout << " Twisting hessian:   " << Kt.norm() << endl;
      }
      
      K = Ks + Kb + Kt;
    }

    // compute collision forces and stiffnesses
    const int rank = R.size();
    SPARSE_MATRIX collisionC(rank, rank);
    if (_collisionsEnabled)
      computeCollisionResponse(R,K,collisionC);

    // compute the LHS of the residual:
    TIMER rhsTimer("Forming RHS");
    _acceleration = -_alpha[0] * (_position - _positionOld) + _alpha[1] * _velocityOld;
    _temp = _M * _acceleration;

    // compute the RHS of the residual:
    // the -_velocityOld is here to match the Baraff-Witkin formulation, and comes
    // out of their variational form for damping. In the direct Rayleigh damping formulation,
    // this extra term would not appear.
    _velocity = _alpha[3] * (_position - _positionOld) - _velocityOld;
    _b = C * _velocity;

    // assemble full residual: LHS + RHS + R - F
    _b += _temp + R + _externalForces;
    rhsTimer.stop();

    // store so we can detect contact breaking later
    unfiltered = _b;

    TIMER lhsTimer("Forming Initial LHS");
    // assemble system matrix A
    //_A = _alpha[0] * _M - _alpha[3] * (C + collisionC) - K;
    _A = (_alpha[0] * _M - _alpha[3] * (C + collisionC) - K).pruned();
    lhsTimer.stop();

    // in [TJM15], this is c = b - Az (page 8, top of column 2)
    TIMER rhsProjectionTimer("RHS PPCG projection");
    VECTOR c = _b - _A * z;
   
    // just for clarity, go ahead and form the RHS and LHS explicitly
    VECTOR rhs = _S * c;
    rhsProjectionTimer.stop();

#if 1
    SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
#else
    /*
    //SPARSE_MATRIX AN = (_A * _N).pruned();
    SPARSE_MATRIX AN = (_A * _N);
    SPARSE_MATRIX ANT = AN.transpose();
  
    //SPARSE_MATRIX leftRight = (_N * AN).pruned();
    SPARSE_MATRIX leftRight = (_N * AN);

    _A += -(AN + ANT) + leftRight + _N;   // final add takes the longest
    const SPARSE_MATRIX& LHS = _A;
    */
    //SPARSE_MATRIX LHS = (_I - _N) * _A * (_I - _N) + _N;
    /*
    SPARSE_MATRIX LHS = (_I - _N) * _A -
                        (_I - _N) * _A * _N.transpose() 
                        + _N.transpose();
                        */

    /*
    SPARSE_MATRIX AN = (_A * _N).pruned();
    SPARSE_MATRIX LHS = -_N * _A -
                        _I * AN +
                        _N * AN;
    LHS = LHS + _A;
    LHS = LHS + _N;
    */
    /*
    SPARSE_MATRIX Nsym = _N - SPARSE_MATRIX(_N.transpose());
    cout << " N sym? " << Nsym.norm() << endl;
    SPARSE_MATRIX Asym = _A - SPARSE_MATRIX(_A.transpose());
    cout << " A sym? " << Asym.norm() << endl;
    */

    TIMER lhsProjectionTimer("LHS PPCG projection");
    SPARSE_MATRIX AN = (_A * _N).pruned();
    // seems like it should just be AN.transpose(), but that bombs
    SPARSE_MATRIX NA = (_N * _A).pruned();
    SPARSE_MATRIX& LHS = _A;
    lhsProjectionTimer.stop();

    TIMER leftRightTimer("LHS left-right projection");
    SPARSE_MATRIX leftRight = (_N * AN).pruned(); // < 1% of running time
    leftRightTimer.stop();

    TIMER lhsFinalAdd("LHS final add");
    LHS = LHS - NA - AN + _N + leftRight;// ~5% of running time
    lhsFinalAdd.stop();
#if 0
    SPARSE_MATRIX ground = _S * _A * _S + _IminusS;
    SPARSE_MATRIX diff = LHS - ground;

    cout << " matrix diff: " << diff.norm() << endl;
    cout << " ground: " << ground.norm() << endl;
    cout << " LHS:    " << LHS.norm() << endl;
    if (diff.norm() > 1e-4)
    {
      for (int k=0; k<diff.outerSize(); ++k)
        for (SPARSE_MATRIX::InnerIterator it(diff,k); it; ++it)
        {
          if (fabs(it.value()) > 1e-4)
          {
            cout << " diff = " << it.value() << " at (" << it.row()<< ", " << it.col() << ")" << endl;   // row index
            //const int row = it.row();
            //const int col = it.col();
            //cout << " N:     " << _N.coeffRef(row, col) << endl;
            //cout << " I - S: " << 1.0 - _N.coeffRef(row, col) << endl;
            //cout << " S:     " << _S.coeffRef(row, col) << endl;
          }
        }
    }
    assert(diff.norm() < 1e-4);
#endif
#endif

    maxR = rhs.squaredNorm();
    if (verbose)
    {
      cout << "done. " << endl;
      cout << " Newton residual: " << maxR << endl;
    }

    VECTOR y = solveSystem(LHS, rhs);

    // aliasing _solution to \Delta x just to make clear what we're doing here
    VECTOR& xDelta = _solution;
    xDelta = y + z;
  
    // update positions
    _position += xDelta;

    // when checking against normals, unfiltered should be negated for Newmark
    //unfiltered *= -1.0;
    //const bool constraintsChanged = findSeparatingSurfaceConstraints(unfiltered);
    const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);

    // see if any of the constraints changed. Used to be that this was outside the Newton loop
    // because the behavior was too oscillatory, but causes too many penetrations to slip
    // through when the Poisson's ratio gets high
    if (constraintsChanged)
    {
      deleteSurfaceConstraints(verbose);
      updateSurfaceConstraints();
      buildConstraintMatrixFaster();
      updateConstraintTargets();
    }
    // update the targets, but the constraint matrix should not have changed.
    else
    {
      updateSurfaceConstraints();
      updateConstraintTargets();
    }
    // update node positions
    _strandMesh.setDisplacement(_position);
   
    step++;
  }
	
  // update velocity
  _velocity = _alpha[3] * (_position - _positionOld);
  //_velocity = (1.0 / _dt) * (_position - _positionOld);

	// update acceleration
  _acceleration = _alpha[0] * (_position - _positionOld) - _alpha[1] * _velocityOld;
  //_acceleration = (1.0 / _dt) * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;
  _acceleration = _S * _acceleration;

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;

  return true;
}

} // HOBAK
} // TIMESTEPPER
