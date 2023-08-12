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
#include "VOLUME_BDF_1.h"
#include "TIMER.h"

#define DEBUG_KINEMATICS 0
using namespace std;

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
VOLUME_BDF_1::VOLUME_BDF_1(TET_STRAND_MESH& strandMesh, 
                           STRAND::STRETCHING& stretching,
                           ENERGY_12D& eeEnergy) :
  VOLUME_TIMESTEPPER(strandMesh, stretching, eeEnergy)
{
  _accelerationOld.resize(_DOFs);
  _acceleration.resize(_DOFs);
  _accelerationOld.setZero();
  _acceleration.setZero();

  _maxNewtonIterations = 3;

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
void VOLUME_BDF_1::setDt(const REAL dt)
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
bool VOLUME_BDF_1::solveDynamics(const bool verbose)
{
  const MATRIX3 R = MATRIX3::Identity();
  const VECTOR3 t = VECTOR3::Zero();

  cout << " Rayleigh: " << _rayleighAlpha << " " << _rayleighBeta << endl;

  return solveDynamicsWithRotation(verbose, R, t);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool VOLUME_BDF_1::solveDynamicsWithRotation(const bool verbose, const MATRIX3& rotation, const VECTOR3& translation)
{
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

  // DEBUG: does it work better here? 
  applyRigidRotation(rotation, translation);

  while (step < _maxNewtonIterations && maxR > eps)
  {
    // build new constraints and see if we should break any
    findNewSurfaceConstraints(verbose);
    VOLUME_TIMESTEPPER::buildConstraintMatrixFaster();
    //VOLUME_TIMESTEPPER::buildConstraintMatrix();

    TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
    if (tetStrandMesh == NULL)
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      std::cout << " Wrong STRAND_MESH type passed to VOLUME_TIMESTEPPER!!!" << std::endl;
      assert(tetStrandMesh);
    }
    tetStrandMesh->computeFs();
    tetStrandMesh->computeSVDs();

    // do collision detection, including spatial data structure updates 
    if (_collisionsEnabled)
      _strandMesh.computeEdgeEdgeCollisions(_edgeEdgeEnergy.eps(), verbose);

    // "z is a vector of the desired values for the constrained variables",
    // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
    // _constraintTargets did not project off the null directions
    updateConstraintTargets();
    z =_IminusS * _constraintTargets;

    // get the internal forces
    TIMER forceTimer("Internal forces");
    const VECTOR stretchingForce = tetStrandMesh->computeStretchingForces();
    const VECTOR bendingForce = tetStrandMesh->computeBendingForces();
    const VECTOR twistingForce = tetStrandMesh->computeTwistingForces();

    if (verbose)
    {
      cout << " Bending force :   " << bendingForce.norm() << endl;
      cout << " Stretching force: " << stretchingForce.norm() << endl;
      cout << " Twisting force:   " << twistingForce.norm() << endl;
      cout << " Solution norm: " << _solution.norm() << endl;

      //cout << endl;
      //cout << " Stretching energy: " << tetStrandMesh->computeStretchingEnergy() << endl;
      //cout << " Bending energy: " << tetStrandMesh->computeBendingEnergy() << endl;
      //cout << " Twisting energy: " << tetStrandMesh->computeTwistingEnergy() << endl;
      cout << " Forming linear system ... " << flush;
    }

    // making non-const so we can add collision forces
    VECTOR R = bendingForce + stretchingForce + twistingForce;
    forceTimer.stop();

    // get the reduced stiffness matrix
    SPARSE_MATRIX K = tetStrandMesh->computeHyperelasticClampedHessianFast();
    //SPARSE_MATRIX K = tetStrandMesh->computeHyperelasticClampedHessian();

    // DEBUG: need better damping here?
    //C = _rayleighAlpha * _M + _rayleighBeta * K;

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

    // assemble system matrix A
    TIMER lhsTimer("Forming LHS");
    _A = _alpha[0] * _M - _alpha[3] * (C + collisionC) - K;

    // in [TJM15], this is c = b - Az (page 8, top of column 2)
    VECTOR c = _b - _A * z;
   
    // just for clarity, go ahead and form the RHS and LHS explicitly
    VECTOR rhs = _S * c;

#define VERIFY_FILTER_BDF 0
#if VERIFY_FILTER_BDF
    SPARSE_MATRIX AN = (_A * _N).pruned();
    SPARSE_MATRIX ANT = AN.transpose();
  
    SPARSE_MATRIX leftRight = (_N * AN).pruned();

    SPARSE_MATRIX ground = _A;
    ground += -(AN + ANT) + leftRight + _N;   // final add takes the longest
#endif
    SPARSE_MATRIX& LHS = filteredSystem();
    lhsTimer.stop();
#if VERIFY_FILTER_BDF
    const SPARSE_MATRIX diff = ground - LHS;
    const REAL diffNorm = diff.norm() / LHS.nonZeros();
    cout << " Filtered system diff: " << diffNorm << endl;
    if (diffNorm > 1e-8)
    {
      cout << " ground: " << endl << clampSmalls(MATRIX(ground)) << endl;
      cout << " LHS: " << endl << clampSmalls(MATRIX(LHS)) << endl;
      cout << " diff: " << endl << clampSmalls(MATRIX(diff)) << endl;
      exit(0);
    }
    else
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " FILTERED SYSTEM MATCHES " << endl;
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    }
    assert(diffNorm < 1e-8);
#endif

#if DEBUG_KINEMATICS
    cout << endl;
    cout << " A eigenvalues: " << endl << eigenvalues(_A).transpose() << endl;
    cout << " K eigenvalues: " << endl << -1.0 * eigenvalues(K).transpose() << endl;
    cout << " LHS eigenvalues: " << endl << eigenvalues(LHS).transpose() << endl;
#endif

    maxR = rhs.squaredNorm();
    if (verbose)
    {
      cout << "done. " << endl;
      cout << " Newton residual: " << maxR << endl;

      if (rhs.hasNaN())
      {
        cout << " RHS contains a NaN!" << endl;
        cout << " c hasNan: " << c.hasNaN() << endl;
        cout << " z hasNan: " << z.hasNaN() << endl;
        cout << " b hasNan: " << _b.hasNaN() << endl;
        cout << " R hasNan: " << R.hasNaN() << endl;
        cout << " velocity hasNan: " << _velocity.hasNaN() << endl;
        cout << " C hasNan: " << MATRIX(C).hasNaN() << endl;
        exit(0);
      }
    }

    VECTOR y = solveSystem(LHS, rhs);
    /*
    cout << " position:      " << _position.norm() << endl;
    cout << " velocity:      " << _velocity.norm() << endl;
    cout << " acceleration:  " << _acceleration.norm() << endl;
    cout << " A norm:        " << MATRIX(_A).squaredNorm() << endl;
    cout << " LHS norm:      " << MATRIX(LHS).squaredNorm() << endl;
    cout << " b norm:        " << _b.norm() << endl;
    cout << " z norm:        " << z.norm() << endl;
    cout << " rhs norm:      " << rhs.norm() << endl;
    cout << " Solution norm: " << y.norm() << endl;
    cout << " LHS: " << endl << MATRIX(LHS) << endl;
    cout << " rhs: " << endl << rhs << endl;
    cout << " sanity: " << (LHS * y - rhs).norm() << endl;
    cout << " solution: " << endl << y << endl;
    */

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
      VOLUME_TIMESTEPPER::buildConstraintMatrixFaster();
      //VOLUME_TIMESTEPPER::buildConstraintMatrix();
      updateConstraintTargets();
    }
    // update the targets, but the constraint matrix should not have changed.
    else
    {
      updateSurfaceConstraints();
      updateConstraintTargets();
    }
    // update node positions
    tetStrandMesh->setDisplacement(_position);
   
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
