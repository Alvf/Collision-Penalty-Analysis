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
#include "NET_TIMESTEPPER.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Geometry/STRAND_NET_MESH.h"
#include "TIMER.h"
#include <float.h>
#include <iostream>
#include "PCG.h"
#include "DIAGONAL.h"
#include "STRAND_DIAGONAL.h"
#include "Hyperelastic/Volume/ANISOTROPIC_ARAP.h"
#include "Hyperelastic/Volume/ANISOTROPIC_DIRICHLET.h"

using namespace std;

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
NET_TIMESTEPPER::NET_TIMESTEPPER(STRAND_NET_MESH& strandMesh, 
                                       STRAND::STRETCHING& stretching) :
  TIMESTEPPER(strandMesh, stretching)
{
  STRAND_NET_MESH* tetStrandMesh = dynamic_cast<STRAND_NET_MESH*>(&_strandMesh);
  if (tetStrandMesh == NULL)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " Wrong STRAND_MESH type passed to NET_TIMESTEPPER!!!" << std::endl;
    assert(tetStrandMesh);
  }
  const VECTOR displacement = tetStrandMesh->getDisplacement();
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " displacement size: " << displacement.size() << endl;
  _position = displacement;
  _positionOld = _position;

  _velocity.resize(_DOFs);
  _velocityOld.resize(_DOFs);

  _velocity.setZero();
  _velocityOld.setZero();
  _velocityDelta.setZero();
  _temp.setZero();
  _solution.setZero();
  _positionOld.setZero();
}

NET_TIMESTEPPER::~NET_TIMESTEPPER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the damping matrix based on the rest pose stiffness
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX NET_TIMESTEPPER::buildRayleighDampingMatrix()
{
  TIMER functionTimer(__FUNCTION__);
  // back up current state
  _temp = _strandMesh.getDisplacement();

  // set to zero displacement
  VECTOR zero(_DOFs);
  zero.setZero();
  _strandMesh.setDisplacement(zero);

  // get stiffness matrix at that state
  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  tetStrandMesh->computeFs();
  tetStrandMesh->computeSVDs();
  SPARSE_MATRIX K = tetStrandMesh->computeHyperelasticClampedHessianFast();

  // restore old state
  _strandMesh.setDisplacement(_temp);
  
  // build out the Rayleigh damping
  SPARSE_MATRIX C = _rayleighAlpha * _M;
  C += _rayleighBeta * K;
  return C;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool NET_TIMESTEPPER::solveDynamics(const bool verbose)
{
  const MATRIX3 R = MATRIX3::Identity();
  const VECTOR3 t = VECTOR3::Zero();
  return solveDynamicsWithRotation(verbose, R, t);
}
#if 0
{
  TIMER functionTimer(__FUNCTION__);
  const bool veryVerbose = false;
  //const bool veryVerbose = true;

  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " STRAND::NET_TIMESTEPPER RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
    cout << " E = " << _strandMesh.E() << endl;
    if (_collisionsEnabled)
      cout << " Collisions are: ON " << endl;
    else
      cout << " Collisions are: OFF " << endl;
  }

  //TIMER prologueTimer("Prologue");

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // build new constraints and see if we should break any
  findNewSurfaceConstraints(verbose);
  //buildConstraintMatrix();
  if (verbose)
  {
    cout << " Kinematic constraints: " << _kinematicConstraints.size() << endl;
    cout << " Plane constraints:     " << _planeConstraints.size() << endl;
  }
  buildConstraintMatrixFaster();

  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  if (tetStrandMesh == NULL)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " Wrong STRAND_MESH type passed to NET_TIMESTEPPER!!!" << std::endl;
    assert(tetStrandMesh);
  }
  tetStrandMesh->computeFs();
  tetStrandMesh->computeSVDs();

  //_strandMesh.setDisplacement(_position);

  // do collision detection, including spatial data structure updates 
  if (_collisionsEnabled)
    _strandMesh.computeEdgeEdgeCollisions(verbose);

  // TODO: add plane constraints so that this starts firing
  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  updateConstraintTargets();
  VECTOR z =_IminusS * _constraintTargets;
  //prologueTimer.stop();

#if 0
  VECTOR R = tetStrandMesh->computeHyperelasticForces();
#else
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
  }

  // making non-const so we can add collision forces
  VECTOR R = bendingForce + stretchingForce + twistingForce;
  forceTimer.stop();
#endif

  SPARSE_MATRIX K = tetStrandMesh->computeHyperelasticClampedHessianFast();

  // collision damping only appears on LHS
  const int rank = R.size();

  SPARSE_MATRIX collisionC(rank,rank);
  if (_collisionsEnabled)
  {
    // compute collision forces and stiffnesses
    computeCollisionResponse(R,K,collisionC, true);
  }

  if (veryVerbose)
  {
    cout << " M size: " << _M.rows() << ", " << _M.cols() << endl;
    cout << " v size: " << _velocity.size() << endl;
    cout << " R size: " << R.size() << endl;
    cout << " f size: " << _externalForces.size() << endl;
  }

  // compute the RHS of the residual:
  //TIMER rhsTimer("Forming Initial RHS");
  const REAL invDt = 1.0 / _dt;
  const REAL invDt2 = invDt * invDt;
  _b = (invDt * _M - C) * _velocity + R + _externalForces;

  // collisionC does *not* appear here, since the damping only depends on
  // v^{t+1}, which only depends on \Delta x, which is the variable
  // we are solving for
  //_b = (invDt * _M - C - collisionC) * _velocity + R + _externalForces;
  //rhsTimer.stop();

  // assemble system matrix A
  TIMER lhsTimer("Forming Initial LHS");
  // TODO: implement collisions
  //_A = _M * invDt2 - collisionC * invDt - K; // w/ Rayleigh
  _A = _M * invDt2 - (C + collisionC) * invDt - K; // w/ Rayleigh
  lhsTimer.stop();

  TIMER rhsProjectionTimer("RHS PPCG projection");
  // in [TJM15], this is c = b - Az (page 8, top of column 2)
  // TODO: add plane constraints
  VECTOR c = _b - _A * z;

  // just for clarity, go ahead and form the RHS and LHS explicitly
  //
  // Since _S is sparse, this multipy could be accelerated significantly, 
  // but leaving it as it is for now
  //VECTOR rhs = _S * c;
  Eigen::VectorXd rhs = _S * c;
  rhsProjectionTimer.stop();

  TIMER lhsProjectionTimer("LHS PPCG projection");

  // DEBUG
#define VERIFY_FILTER 0
#if VERIFY_FILTER
  SPARSE_MATRIX AN = (_A * _N).pruned();
  SPARSE_MATRIX ANT = AN.transpose();
  
  SPARSE_MATRIX leftRight = (_N * AN).pruned();

  SPARSE_MATRIX ground = _A;
  ground += -(AN + ANT) + leftRight + _N;   // final add takes the longest
#endif
  
  const SPARSE_MATRIX& LHS = filteredSystem();
#if VERIFY_FILTER
  const SPARSE_MATRIX diff = ground - LHS;
  const REAL diffNorm = diff.norm() / LHS.nonZeros();
  cout << " Filtered system diff: " << diffNorm << endl;
  if (diffNorm > 1e-8)
  {
    cout << " ground: " << endl << MATRIX(ground) << endl;
    cout << " LHS: " << endl << MATRIX(LHS) << endl;
    cout << " diff: " << endl << MATRIX(diff) << endl;
  }
  assert(diffNorm < 1e-8);
#endif

  lhsProjectionTimer.stop();

  //if (verbose)
  if (veryVerbose)
  {
    cout << " b norm: " << _b.norm() << endl;
    cout << " r norm: " << c.norm() << endl;
    cout << " z norm: " << z.norm() << endl;
    cout << " targets norm: " << _constraintTargets.norm() << endl;
    cout << " rhs norm: " << rhs.norm() << endl;
    cout << " LHS norm: " << LHS.norm() << endl;
    cout << " A norm:   " << _A.norm() << endl;
    cout << " K norm:   " << K.norm() << endl;
    cout << " R norm:   " << R.norm() << endl;
    cout << " C norm:   " << C.norm() << endl;
  }

  // solve system using whatever solver is activated right now
  VECTOR y = solveSystem(LHS, rhs);
 
  TIMER postTimer("Epilogue");

  // aliasing _solution to \Delta x just to make clear what we're doing here
  VECTOR& xDelta = _solution;
  // TODO: add plane constraints
  xDelta = y + z;
  //xDelta = y;

  if (veryVerbose)
  {
    cout << " position: " << _position.size() << endl;
    cout << " x delta: " << xDelta.size() << endl;
  }

  // update positions
  _position += xDelta;

  if (verbose)
    cout << " xDelta norm: " << xDelta.norm() << endl;

  // TODO: add plane constraints
  // when checking against normals, unfiltered should be negated for Newmark
  const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);
  //const bool constraintsChanged = findSeparatingSurfaceConstraints(_A * xDelta - _b);

  // see if any of the constraints changed. Used to be that this was outside the Newton loop
  // because the behavior was too oscillatory, but causes too many penetrations to slip
  // through when the Poisson's ratio gets high
  if (constraintsChanged)
  {
    deleteSurfaceConstraints(verbose);
    updateSurfaceConstraints();
    //buildConstraintMatrix();
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

  // update velocity
  _velocity = invDt * (_position - _positionOld);

	// update acceleration
  //_acceleration = invDt * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;

  postTimer.stop();

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;
  return true;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool NET_TIMESTEPPER::solveDynamicsWithRotation(const bool verbose, const MATRIX3& rotation, const VECTOR3& translation)
{
  TIMER functionTimer(__FUNCTION__);
  const bool veryVerbose = false;
  //const bool veryVerbose = true;

  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " STRAND::NET_TIMESTEPPER RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
    cout << " E = " << _strandMesh.E() << endl;
    if (_collisionsEnabled)
      cout << " Collisions are: ON " << endl;
    else
      cout << " Collisions are: OFF " << endl;
  }

  //TIMER prologueTimer("Prologue");

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;

  // DEBUG: does it work better here? 
  applyRigidRotation(rotation, translation);

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // build new constraints and see if we should break any
  findNewSurfaceConstraints(verbose);
  if (verbose)
  {
    cout << " Kinematic constraints: " << _kinematicConstraints.size() << endl;
    cout << " Plane constraints:     " << _planeConstraints.size() << endl;
  }
  //buildConstraintMatrix();
  buildConstraintMatrixFaster();

  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  if (tetStrandMesh == NULL)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " Wrong STRAND_MESH type passed to NET_TIMESTEPPER!!!" << std::endl;
    assert(tetStrandMesh);
  }
  tetStrandMesh->computeFs();
  tetStrandMesh->computeSVDs();

  //_strandMesh.setDisplacement(_position);

  // do collision detection, including spatial data structure updates 
  if (_collisionsEnabled)
    _strandMesh.computeEdgeEdgeCollisions(verbose);

  // TODO: add plane constraints so that this starts firing
  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  updateConstraintTargets();
  VECTOR z =_IminusS * _constraintTargets;
  //prologueTimer.stop();

#if 0
  VECTOR R = tetStrandMesh->computeHyperelasticForces();
#else
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
  }

  // making non-const so we can add collision forces
  VECTOR R = bendingForce + stretchingForce + twistingForce;
  forceTimer.stop();
#endif

  SPARSE_MATRIX K = tetStrandMesh->computeHyperelasticClampedHessianFast();

  // collision damping only appears on LHS
  const int rank = R.size();

  SPARSE_MATRIX collisionC(rank,rank);
  if (_collisionsEnabled)
  {
    // compute collision forces and stiffnesses
    computeCollisionResponse(R,K,collisionC, true);
  }

  if (veryVerbose)
  {
    cout << " M size: " << _M.rows() << ", " << _M.cols() << endl;
    cout << " v size: " << _velocity.size() << endl;
    cout << " R size: " << R.size() << endl;
    cout << " f size: " << _externalForces.size() << endl;
  }

  // compute the RHS of the residual:
  //TIMER rhsTimer("Forming Initial RHS");
  const REAL invDt = 1.0 / _dt;
  const REAL invDt2 = invDt * invDt;
  _b = (invDt * _M - C) * _velocity + R + _externalForces;

  // collisionC does *not* appear here, since the damping only depends on
  // v^{t+1}, which only depends on \Delta x, which is the variable
  // we are solving for
  //_b = (invDt * _M - C - collisionC) * _velocity + R + _externalForces;
  //rhsTimer.stop();

  // assemble system matrix A
  TIMER lhsTimer("Forming Initial LHS");
  // TODO: implement collisions
  //_A = _M * invDt2 - collisionC * invDt - K; // w/ Rayleigh
  _A = _M * invDt2 - (C + collisionC) * invDt - K; // w/ Rayleigh
  lhsTimer.stop();

  TIMER rhsProjectionTimer("RHS PPCG projection");
  // in [TJM15], this is c = b - Az (page 8, top of column 2)
  // TODO: add plane constraints
  VECTOR c = _b - _A * z;

  // just for clarity, go ahead and form the RHS and LHS explicitly
  //
  // Since _S is sparse, this multipy could be accelerated significantly, 
  // but leaving it as it is for now
  //VECTOR rhs = _S * c;
  Eigen::VectorXd rhs = _S * c;
  rhsProjectionTimer.stop();

  TIMER lhsProjectionTimer("LHS PPCG projection");

  // DEBUG
#define VERIFY_FILTER 0
#if VERIFY_FILTER
  SPARSE_MATRIX AN = (_A * _N).pruned();
  SPARSE_MATRIX ANT = AN.transpose();
  
  SPARSE_MATRIX leftRight = (_N * AN).pruned();

  SPARSE_MATRIX ground = _A;
  ground += -(AN + ANT) + leftRight + _N;   // final add takes the longest
#endif
  
  const SPARSE_MATRIX& LHS = filteredSystem();
#if VERIFY_FILTER
  const SPARSE_MATRIX diff = ground - LHS;
  const REAL diffNorm = diff.norm() / LHS.nonZeros();
  cout << " Filtered system diff: " << diffNorm << endl;
  if (diffNorm > 1e-8)
  {
    cout << " ground: " << endl << MATRIX(ground) << endl;
    cout << " LHS: " << endl << MATRIX(LHS) << endl;
    cout << " diff: " << endl << MATRIX(diff) << endl;
  }
  assert(diffNorm < 1e-8);
#endif

  lhsProjectionTimer.stop();

  if (verbose)
  // if (veryVerbose)
  {
    cout << " b norm: " << _b.norm() << endl;
    cout << " r norm: " << c.norm() << endl;
    cout << " z norm: " << z.norm() << endl;
    cout << " targets norm: " << _constraintTargets.norm() << endl;
    cout << " rhs norm: " << rhs.norm() << endl;
    cout << " LHS norm: " << LHS.norm() << endl;
    cout << " A norm:   " << _A.norm() << endl;
    cout << " K norm:   " << K.norm() << endl;
    cout << " R norm:   " << R.norm() << endl;
    cout << " C norm:   " << C.norm() << endl;
  }

  // solve system using whatever solver is activated right now
  VECTOR y = solveSystem(LHS, rhs);
 
  TIMER postTimer("Epilogue");

  // aliasing _solution to \Delta x just to make clear what we're doing here
  VECTOR& xDelta = _solution;
  // TODO: add plane constraints
  xDelta = y + z;
  //xDelta = y;

  if (veryVerbose)
  {
    cout << " position: " << _position.size() << endl;
    cout << " x delta: " << xDelta.size() << endl;
  }

  // update positions
  _position += xDelta;

  if (verbose)
    cout << " xDelta norm: " << xDelta.norm() << endl;

  // TODO: add plane constraints
  // when checking against normals, unfiltered should be negated for Newmark
  const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);
  //const bool constraintsChanged = findSeparatingSurfaceConstraints(_A * xDelta - _b);

  // see if any of the constraints changed. Used to be that this was outside the Newton loop
  // because the behavior was too oscillatory, but causes too many penetrations to slip
  // through when the Poisson's ratio gets high
  if (constraintsChanged)
  {
    deleteSurfaceConstraints(verbose);
    updateSurfaceConstraints();
    //buildConstraintMatrix();
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

  // update velocity
  _velocity = invDt * (_position - _positionOld);

	// update acceleration
  //_acceleration = invDt * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;

  postTimer.stop();

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void NET_TIMESTEPPER::stiffenKinematicStrands()
{
  // find all the constrained vertices
  map<int, bool> vertexConstrained;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    vertexConstrained[_kinematicConstraints[x].vertexID] = true;

  // look at all the edges
  map<int, bool> edgeConstrained;
  const vector<VECTOR2I>& edgeIndices = _strandMesh.edgeIndices();
  for (unsigned int x = 0; x < edgeIndices.size(); x++)
  {
    bool first  = false;
    bool second = false;

    // is one end constrained while the other isn't?
    if (vertexConstrained.find(edgeIndices[x][0]) != vertexConstrained.end())
      first = true;
    if (vertexConstrained.find(edgeIndices[x][1]) != vertexConstrained.end())
      second = true;

    // if neither is constrained, move on
    //if (!first && !second)
    // XOR: only care if one is constrained and the other is not
    if (!(first ^ second))
      continue;

    // tag it to be stiffened
    edgeConstrained[x] = true;
  }

  // look at all the tets
  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  const vector<VECTOR4I>& tets = tetStrandMesh->tets();
  const vector<VECTOR3I>& tetEdges = tetStrandMesh->tetEdges();
  vector<STRAND::COMPOSITE>& materials = tetStrandMesh->materials();
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    // should one of its edges be stiffened?
    const VECTOR3I& edges = tetEdges[x];

    bool stiffen = false;
    for (int y = 0; y < 3; y++)
      if (edgeConstrained.find(edges[y]) != edgeConstrained.end())
        stiffen = true;

    if (!stiffen) continue;
    //cout << " Reinforcing tet " << x << endl;

    // stiffen the twisting energy
    vector<VOLUME::HYPERELASTIC*>& volumeEnergies = materials[x].volumeEnergies();
    for (unsigned int y = 0; y < volumeEnergies.size(); y++)
    {
      using namespace VOLUME;
      ANISOTROPIC_ARAP* arap = dynamic_cast<ANISOTROPIC_ARAP*>(volumeEnergies[y]);
      if (arap != NULL)
        arap->mu() *= 100.0;
      /*
      ANISOTROPIC_DIRICHLET* dirichlet = dynamic_cast<ANISOTROPIC_DIRICHLET*>(volumeEnergies[y]);
      if (dirichlet != NULL)
        dirichlet->mu() *= 100.0;
        */
    }

    /*
    // stiffen the bending energy
    std::array<ISOTROPIC_BENDING*,2>& bendingEnergies = materials[x].bendingEnergies();
    if (bendingEnergies[0] != NULL)
      bendingEnergies[0]->mu() *= 2.0;
    if (bendingEnergies[1] != NULL)
      bendingEnergies[1]->mu() *= 2.0;
      */
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX NET_TIMESTEPPER::buildPlaneConstraintsOnly()
{
  TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsP;

  vector<bool> diagonalSeen(_DOFs);
  for (int x = 0; x < _DOFs; x++)
    diagonalSeen[x] = false;

  vector<bool> isKinematic(_strandMesh.totalVertices());
  for (int x = 0; x < _strandMesh.totalVertices(); x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

  // build the plane constraints for the LHS
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // if this one is tagged for deletion, ignore it
    if (constraint.isSeparating) continue;

    // get the normal direction
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const VECTOR3& localNormal = _planeConstraints[x].localNormal;
    const VECTOR3 normal = shape->localNormalToWorld(localNormal).normalized();

    // build the filter matrix
    const MATRIX3 Nblock = normal * normal.transpose();
    const MATRIX3 Sblock = MATRIX3::Identity() - Nblock;
    const int vertexID = constraint.vertexID;

    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    const int index = 3 * vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        if (index + i == index + j) 
          diagonalSeen[index + i] = true;
        TRIPLET tripletP(index + i, index + j, Sblock(i,j));
        tripletsP.push_back(tripletP);
      }
  }

  // if the diagonal was never set, set it to one
  for (int x = 0; x < _DOFs; x++)
  {
    if (diagonalSeen[x]) continue;
    TRIPLET tripletP(x, x, 1);
    tripletsP.push_back(tripletP);
  }

  SPARSE_MATRIX P(_DOFs, _DOFs);
  P.setFromTriplets(tripletsP.begin(), tripletsP.end());

  return P;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX NET_TIMESTEPPER::buildPlaneConstraintsNoIdentity()
{
  // timing is negligible
  //TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsP;

  vector<bool> isKinematic(_strandMesh.totalVertices());
  for (int x = 0; x < _strandMesh.totalVertices(); x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

  // build the plane constraints for the LHS
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // if this one is tagged for deletion, ignore it
    if (constraint.isSeparating) continue;

    // get the normal direction
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const VECTOR3& localNormal = _planeConstraints[x].localNormal;
    const VECTOR3 normal = shape->localNormalToWorld(localNormal).normalized();

    // build the filter matrix
    const MATRIX3 Nblock = normal * normal.transpose();
    const MATRIX3 Sblock = MATRIX3::Identity() - Nblock;
    const int vertexID = constraint.vertexID;

    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    const int index = 3 * vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        TRIPLET tripletP(index + i, index + j, Sblock(i,j));
        tripletsP.push_back(tripletP);
      }
  }

  SPARSE_MATRIX P(_DOFs, _DOFs);
  P.setFromTriplets(tripletsP.begin(), tripletsP.end());

  return P;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void NET_TIMESTEPPER::buildConstraintMatrix()
{
  TIMER functionTimer(__FUNCTION__);
  SPARSE_MATRIX I(_DOFs, _DOFs);
  I.setIdentity();
  _S = I;

  vector<bool> isKinematic(_strandMesh.totalVertices());
  for (int x = 0; x < _strandMesh.totalVertices(); x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

#if 0
  // output all the constrained vertexIDs
  const bool debug = false;
  if (debug)
  {
    const vector<VECTOR3>& vertices = _strandMesh.vertices();
    cout << " Plane constrained vertices: " << flush;
    for (unsigned int x = 0; x < _planeConstraints.size(); x++)
    {
      const PLANE_CONSTRAINT& constraint = _planeConstraints[x];
      const int vertexID = constraint.vertexID;
      cout << " Vertex " << vertexID << ": " << vertices[vertexID].transpose() << endl;
    }
    cout << endl;
    
    cout << " Pinned vertices: " << flush;
    for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    {
      const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];
      const int vertexID = constraint.vertexID;
      cout << " Vertex " << vertexID << ": " << vertices[vertexID].transpose() << endl;
    }
    cout << endl;
  }
#endif

  // build the plane constraints for the LHS
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // if this one is tagged for deletion, ignore it
    if (constraint.isSeparating) continue;

    // get the normal direction
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const VECTOR3& localNormal = _planeConstraints[x].localNormal;
    const VECTOR3 normal = shape->localNormalToWorld(localNormal).normalized();

    // build the filter matrix
    const MATRIX3 Sblock = MATRIX3::Identity() - normal * normal.transpose();
    const int vertexID = constraint.vertexID;

    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    const int index = 3 * vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        _S.coeffRef(index + i, index + j) = Sblock(i,j);
  }

  // apply the kinematic constraints LAST. These override any prior plane
  // constraints
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // set the filter matrix entries
    const int index = 3 * constraint.vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        _S.coeffRef(index + i, index + j) = 0.0;
  }

 /* 
  unsigned int vertexEnd = 3 * _strandMesh.totalVertices();
  for (unsigned int x = 0; x < _constrainedEdges.size(); x++)
  {
    unsigned int index = vertexEnd + _constrainedEdges[x];
    _S.coeffRef(index,index) = 0.0;
  }
  */

  // store the complement
  _IminusS = I - _S;
  _N = _IminusS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void NET_TIMESTEPPER::buildConstraintMatrixFaster()
{
  TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsS;
  vector<TRIPLET> tripletsN;

  vector<bool> diagonalSeen(_DOFs);
  for (int x = 0; x < _DOFs; x++)
    diagonalSeen[x] = false;

  vector<bool> isKinematic(_strandMesh.totalVertices());
  for (int x = 0; x < _strandMesh.totalVertices(); x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

#if 0
  // output all the constrained vertexIDs
  const bool debug = true;
  if (debug)
  {
    cout << " Plane constrained vertices: " << flush;
    for (unsigned int x = 0; x < _planeConstraints.size(); x++)
    {
      const PLANE_CONSTRAINT& constraint = _planeConstraints[x];
      const int vertexID = constraint.vertexID;
      cout << vertexID << " " << flush;
    }
    cout << endl;
    
    cout << " Pinned vertices: " << flush;
    for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    {
      const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];
      const int vertexID = constraint.vertexID;
      cout << vertexID << " " << flush;
    }
    cout << endl;
  }
#endif

  // build the plane constraints for the LHS
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // if this one is tagged for deletion, ignore it
    if (constraint.isSeparating) continue;

    // get the normal direction
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const VECTOR3& localNormal = _planeConstraints[x].localNormal;
    const VECTOR3 normal = shape->localNormalToWorld(localNormal).normalized();

    // build the filter matrix
    const MATRIX3 Nblock = normal * normal.transpose();
    const MATRIX3 Sblock = MATRIX3::Identity() - Nblock;
    const int vertexID = constraint.vertexID;

    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    const int index = 3 * vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        if (index + i == index + j) 
          diagonalSeen[index + i] = true;
        TRIPLET tripletS(index + i, index + j, Sblock(i,j));
        tripletsS.push_back(tripletS);
        
        TRIPLET tripletN(index + i, index + j, Nblock(i,j));
        tripletsN.push_back(tripletN);
      }
  }

  // apply the kinematic constraints LAST. These override any prior plane
  // constraints
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // set the filter matrix entries
    const int index = 3 * constraint.vertexID;

    // N block is identity
    for (unsigned int i = 0; i < 3; i++)
    {
      TRIPLET tripletN(index + i, index + i, 1);
      tripletsN.push_back(tripletN);
    }

    // S block just gets zeroed out
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        if (index + i == index + j) 
          diagonalSeen[index + i] = true;
        TRIPLET tripletS(index + i, index + j, 0);
        tripletsS.push_back(tripletS);
      }
  }

  // if the diagonal was never set, set it to one
  for (int x = 0; x < _DOFs; x++)
  {
    if (diagonalSeen[x]) continue;
    TRIPLET tripletS(x, x, 1);
    tripletsS.push_back(tripletS);
  }

  _S = SPARSE_MATRIX(_DOFs, _DOFs);
  _S.setFromTriplets(tripletsS.begin(), tripletsS.end());
  
  _N = _I - _S;

  // DEBUG: is this needed?
  //_I = SPARSE_MATRIX(_DOFs, _DOFs);
  //_I.setIdentity();

  // store the complement
  //_IminusS = _I - _S;
  _IminusS = _N;

#if 0
  SPARSE_MATRIX diff = (_N - (_I - _S));
  cout << " I - S Diff: " << diff.norm() << endl;
  if (diff.norm() > 1e-4)
  {
    for (int k=0; k<diff.outerSize(); ++k)
      for (SPARSE_MATRIX::InnerIterator it(diff,k); it; ++it)
      {
        if (fabs(it.value()) > 0)
        {
          cout << " diff = " << it.value() << " at (" << it.row()<< ", " << it.col() << ")" << endl;   // row index
          const int row = it.row();
          const int col = it.col();
          cout << " N:     " << _N.coeffRef(row, col) << endl;
          cout << " I - S: " << 1.0 - _N.coeffRef(row, col) << endl;
          cout << " S:     " << _S.coeffRef(row, col) << endl;
        }
      }
  }
  assert(diff.norm() < 1e-4);
  diff = _S - (_I - _N);
  cout << " S Diff: " << diff.norm() << endl;
  assert(diff.norm() < 1e-4);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// indirect function so we can try out lots of different solvers
///////////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR NET_TIMESTEPPER::solveSystem(const SPARSE_MATRIX& LHS, const VECTOR& rhs)
{
  cout << " Solving system ..." << flush;
  Eigen::VectorXd y = Eigen::VectorXd::Zero(LHS.rows());

  bool verbose = false;
  bool usingAMG = false;
  bool usingCR = false;
  //bool usingEigenCG = false;
  bool usingEigenCG = _disablePreconditioner;

  if (!_pcgEnabled)
  {
    if (verbose)
      cout << " USING CHOLESKY" << endl;
    TIMER choleskyTimer("Cholesky Solve");
    Eigen::SimplicialLDLT<SPARSE_MATRIX> solver;
    solver.compute(LHS);
    y = solver.solve(rhs);
    choleskyTimer.stop();
    cout << " done." << endl;
    return y;
  }

  if (usingCR)
  {
    if (verbose)
      cout << " USING CR" << endl;
    TIMER pcgTimer("CR Solve");

    VECTOR residual = LHS * y - rhs;
    cout << " Residual: " << residual.norm() << endl;
    
    DIAGONAL diagonal(LHS);
    PCG pcgSolver(LHS, diagonal);
    //y = pcgSolver.solveCR(rhs);
    y = pcgSolver.solvePCR(rhs);
    pcgTimer.stop();

    if (verbose)
      printf("  PCR iters: %3i err: %6.4e \n", (int)pcgSolver.iterations(), (float)pcgSolver.error());
  }
  else if (usingEigenCG)
  {
    if (verbose)
      cout << " USING Eigen PCG" << endl;
    TIMER pcgTimer("PCG Solve");
    _cgSolver.compute(LHS);
    y = _cgSolver.solve(rhs);
    pcgTimer.stop();

    if (verbose)
      printf("  Eigen PCG iters: %3i err: %6.4e \n", (int)_cgSolver.iterations(), (float)_cgSolver.error());
  }
  else
  {
    if (verbose)
      cout << " USING PCG" << endl;
    TIMER pcgTimer("PCG Solve");

    VECTOR residual = LHS * y - rhs;
    //cout << " Newton residual: " << residual.norm() << endl;
    
    STRAND_DIAGONAL diagonal(LHS, _strandMesh.globalStrandEnds());
    PCG pcgSolver(LHS, diagonal);
    y = pcgSolver.solveEigenStyle(rhs);
    pcgTimer.stop();

    if (verbose)
      printf("  PCG iters: %3i err: %6.4e \n", (int)pcgSolver.iterations(), (float)pcgSolver.error());
  }
  cout << " done." << endl;
  return y; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// apply PPCG filters to system matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& NET_TIMESTEPPER::filteredSystem()
{
  TIMER functionTimer("PPCG: filteredSystem");
#if 1

  // apply plane contraints
  _P = buildPlaneConstraintsOnly();
  _PnoI = buildPlaneConstraintsNoIdentity();

  // see which kinematic entries are zero
  vector<bool> isKinematic(3 * _strandMesh.totalVertices());
  for (unsigned int x = 0; x < isKinematic.size(); x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    unsigned int index = _kinematicConstraints[x].vertexID;
    isKinematic[3 * index] = true;
    isKinematic[3 * index + 1] = true;
    isKinematic[3 * index + 2] = true;
  }
  
  // see which plane entries are zero
  vector<bool> isPlanar(3 * _strandMesh.totalVertices());
  for (unsigned int x = 0; x < isPlanar.size(); x++)
    isPlanar[x] = false;
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // if this one is tagged for deletion, ignore it
    if (constraint.isSeparating) continue;
    
    // if it's kinematic, move on
    const int vertexID = constraint.vertexID;
    if (isKinematic[3 * vertexID]) continue;

    isPlanar[3 * vertexID] = true;
    isPlanar[3 * vertexID + 1] = true;
    isPlanar[3 * vertexID + 2] = true;
  }

  _BA = _A;

#if 0
  // BLOCK DIAGONAL TERM WAS MISSING FROM PLANE CONSTRAINT
  SPARSE_MATRIX _missing = _PnoI * _A * _PnoI;

  // zero off the plane constraints
  for (int i = 0; i < _A.outerSize(); i++)
  {
    const int k_start = _A.outerIndexPtr()[i];
    const int k_end   = _A.outerIndexPtr()[i+1];

    for (int k = k_start; k < k_end; k++) 
    {
      int j = _A.innerIndexPtr()[k];
      const bool constrained = (isPlanar[i] || isPlanar[j]);
      if (!constrained) continue;

      _A.valuePtr()[k] = 0;

      if (isPlanar[i])
        _BA.valuePtr()[k] = 0;
    }
  }
#else
  // build out just the block diagonal entries
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;

  // zero off the plane constraints
  for (int i = 0; i < _A.outerSize(); i++)
  {
    const int k_start = _A.outerIndexPtr()[i];
    const int k_end   = _A.outerIndexPtr()[i+1];

    for (int k = k_start; k < k_end; k++) 
    {
      int j = _A.innerIndexPtr()[k];
      const bool constrained = (isPlanar[i] || isPlanar[j]);
      if (!constrained) continue;

      // right before we delete the block diagonal from A,
      // store it here so we can multiply against _PnoI
      TRIPLET triplet(i,j, _A.valuePtr()[k]);
      triplets.push_back(triplet);

      _A.valuePtr()[k] = 0;

      // only zero off the row
      if (isPlanar[i])
        _BA.valuePtr()[k] = 0;
    }
  }
  
  // BLOCK DIAGONAL TERM WAS MISSING FROM PLANE CONSTRAINT
  if (_blockDiag.rows() != _DOFs)
    _blockDiag = SPARSE_MATRIX(_DOFs, _DOFs);
  _blockDiag.setFromTriplets(triplets.begin(), triplets.end());
  _missing = _PnoI * _blockDiag * _PnoI;
#endif

  //cout << " wiped: " << endl << MATRIX(_A) << endl;
  TIMER midAdd("PPCG: Mid add"); // 7.9% 5.7% with +=, 7.5% in twist test
  _A += SPARSE_MATRIX(_PnoI * _A * _PnoI).pruned(1e-7) - _P;  // prune doesn't seem to matter, 
                                                              // inline multiply doesn't seem to matter

  /*
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  std::cout << " after PnoI A PnoI - P" << endl;
  cout << " A: " << endl << clampSmalls(MATRIX(_A)) << endl;
  cout << " PnoI A PnoI: " << endl << clampSmalls(_missing) << endl;
  */

  // add identity by hand
//#pragma omp parallel
//#pragma omp for schedule(static)
  for (int i = 0; i < _A.outerSize(); i++)
  {
    const int k_start = _A.outerIndexPtr()[i];
    const int k_end   = _A.outerIndexPtr()[i+1];

    for (int k = k_start; k < k_end; k++) 
    {
      int j = _A.innerIndexPtr()[k];
      if (i == j)
        _A.valuePtr()[k] += 1;
    }
  }
  midAdd.stop();

  _C = (_BA * _PnoI).pruned(1e-7);

  //TIMER addC("PPCG: Add C");  // 6.24% 6.5% 6.12%
  _A += _C + SPARSE_MATRIX(_C.transpose()); // using the temp works well!
  //addC.stop();

  /*
  // everything is back except block diagonals ....
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " BA: " << endl << clampSmalls(MATRIX(_BA)) << endl;
  cout << " PnoI: " << endl << clampSmalls(MATRIX(_PnoI)) << endl;
  cout << " C: " << endl << clampSmalls(MATRIX(_C)) << endl;
  cout << " A: " << endl << clampSmalls(MATRIX(_A)) << endl;
  */

//#pragma omp parallel
//#pragma omp for schedule(static)
  for (int i = 0; i < _A.outerSize(); i++)
  {
    const int k_start = _A.outerIndexPtr()[i];
    const int k_end   = _A.outerIndexPtr()[i+1];

    for (int k = k_start; k < k_end; k++) 
    {
      int j = _A.innerIndexPtr()[k];
      const bool constrained = (isKinematic[i] || isKinematic[j]);
      if (!constrained) continue;

      _A.valuePtr()[k] = (j != i) ? 0.0 : 1.0;
    }
  }
  _A += _missing;
  /*
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " A: " << endl << clampSmalls(MATRIX(_A)) << endl;
  */
  
  return _A;
#else
  SPARSE_MATRIX AN = (_A * _N).pruned();
  SPARSE_MATRIX ANT = AN.transpose();
  
  SPARSE_MATRIX leftRight = (_N * AN).pruned();

  _A += -(AN + ANT) + leftRight + _N;   // final add takes the longest
  return _A;
#endif
}

} // HOBAK
} // STRAND
