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
#include "TIMESTEPPER.h"
#include "TIMER.h"
#include <float.h>
#include <iostream>

using namespace std;

namespace HOBAK {
namespace SHELL {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
TIMESTEPPER::TIMESTEPPER(TRIANGLE_MESH& triangleMesh, 
                         SHELL::STRETCHING& stretching,
                         SHELL::BENDING& bending, 
                         ENERGY_12D& vfEnergy, 
                         ENERGY_12D& eeEnergy) :
  _triangleMesh(triangleMesh), _stretchingEnergy(stretching), _bendingEnergy(bending), _vfEnergy(vfEnergy), _eeEnergy(eeEnergy)
{
  initialize();

  computeRestThetas();
}

TIMESTEPPER::~TIMESTEPPER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::initialize()
{
  _residual = FLT_MAX;
  _seenPCGIterations    = -1;
  _currentTimestep = 0;

  _DOFs = _triangleMesh.DOFs();
  _b.resize(_DOFs);
  _forces.resize(_DOFs);
  _externalForces.resize(_DOFs);
  _constraintTargets.resize(_DOFs);

  _b.setZero();
  _forces.setZero();
  _externalForces.setZero();
  _constraintTargets.setZero();

  int totalVertices = _triangleMesh.vertices().size();
  _inCollision.resize(totalVertices);
  for (int x = 0; x < totalVertices; x++)
    _inCollision[x] = false;

  _position.resize(_DOFs);
  _positionOld.resize(_DOFs);
  _velocity.resize(_DOFs);
  _temp.resize(_DOFs);
  _solution.resize(_DOFs);
  
  _position.setZero();
  _positionOld.setZero();
  _velocity.setZero();
  _temp.setZero();
  _solution.setZero();

  _name = string("UNKNOWN");

  _vertexFaceSelfCollisionsOn = true;
  _edgeEdgeSelfCollisionsOn   = true;
  _collisionStiffness = 1.0;
  _collisionDampingBeta = 0.001;

  _dt = 1.0 / 30.0;

  // build the mass matrix once and for all
  _M = buildMassMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// compute the rest bending for each flap
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeRestThetas()
{
  TIMER functionTimer(__FUNCTION__);
  const vector<VECTOR3>& vertices = _triangleMesh.vertices();
  const vector<VECTOR4I>& flaps = _triangleMesh.flaps();

  vector<REAL>& restThetas = _triangleMesh.restThetas();
  restThetas.clear();
  for (unsigned int i = 0; i < flaps.size(); i++)
  {
    vector<VECTOR3> flapVertices;
    for (unsigned int j = 0; j < 4; j++)
      flapVertices.push_back(vertices[flaps[i][j]]);
    const REAL theta = _bendingEnergy.restAngle(flapVertices);
    restThetas.push_back(theta);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the displacement targets the the Baraff-Witkin-style constraints
// are trying to hit. Assumes that buildConstraintMatrix() has already been called
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateConstraintTargets()
{
  TIMER functionTimer(__FUNCTION__);
  _constraintTargets.setZero();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    // should ignore if we've tagged it for deletion
    if (_planeConstraints[x].isSeparating)
      continue;

    // retrieve collision information
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const int vertexID = constraint.vertexID;
    const int index = 3 * vertexID;

    // compute the target displacement
    const VECTOR3& vertex = _triangleMesh.vertices()[vertexID];
    const VECTOR3& localClosestPoint = _planeConstraints[x].localClosestPoint;
    const VECTOR3& closestPoint = shape->localVertexToWorld(localClosestPoint);

    const VECTOR3& displacement = closestPoint - vertex;
    for (int i = 0; i < 3; i++)
      _constraintTargets[index + i] = displacement[i];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// filter positions to incorporate Baraff-Witkin-style constraints
//
// NOTE: this is different from the quasistatic case, in that it stores the kinematic
// corrections in _position, and not directly the vertices of _tetMesh. When the
// displacement of _tetMesh is pinned during the Newton solve, then the _teMesh
// will see the constraints.
//
// The QUASISTATIC class overrides this one with its own implementation, but the
// dynamics solve all share this version
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::applyKinematicConstraints()
{
  const vector<VECTOR3>& restVertices = _triangleMesh.restVertices();
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // pin the tet mesh position according to the constraint
    const VECTOR3& localPosition = constraint.localPosition;
    VECTOR3 world = constraint.shape->localVertexToWorld(localPosition);

    const int vertexID = constraint.vertexID;
    const VECTOR3 diff = world - restVertices[vertexID];
    _position[3 * vertexID] = diff[0];
    _position[3 * vertexID + 1] = diff[1];
    _position[3 * vertexID + 2] = diff[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the constraint matrix to incorporate Baraff-Witkin-style constraints, but using 
// the "Smoothed aggregation multigrid for cloth simulation" projection matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrix()
{
  SPARSE_MATRIX I(_DOFs, _DOFs);
  I.setIdentity();
  _S = I;

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

  // store the complement
  _IminusS = I - _S;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildMassMatrix()
{
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  const vector<REAL>& areas = _triangleMesh.restOneRingAreas(); 

  // set diagonal the one-ring volumes
  // TODO: take into account density
  for (int x = 0; x < _triangleMesh.totalVertices(); x++)
  {
    const REAL entry = areas[x];
    for (int y = 0; y < 3; y++)
    {
      TRIPLET triplet(3 * x + y, 3 * x + y, entry);
      triplets.push_back(triplet);
    }
  }
  
  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
//
// solves with self collisions, and collision forces are Rayleigh damped
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::solve(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " SHELL::TIMESTEPPER RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

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
  buildConstraintMatrix();

  _triangleMesh.setDisplacement(_position);
  _triangleMesh.computeFs();
  _triangleMesh.computeSVDs();

  // do collision detection, including spatial data structure updates 
  computeCollisionDetection();

  // TODO: add plane constraints so that this starts firing
  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  updateConstraintTargets();
  VECTOR z =_IminusS * _constraintTargets;

  // get the internal forces
  VECTOR R = _triangleMesh.computeStretchingForces(_stretchingEnergy);
  R += _triangleMesh.computeBendingForces(_bendingEnergy);

  // get the reduced stiffness matrix
  SPARSE_MATRIX K = _triangleMesh.computeStretchingClampedHessian(_stretchingEnergy);
  K += _triangleMesh.computeBendingClampedHessian(_bendingEnergy);
#if VERY_VERBOSE
  SPARSE_MATRIX elasticK = -K;
#endif

  // collision damping only appears on LHS
  const int rank = R.size();
  SPARSE_MATRIX collisionC(rank, rank);

  // compute collision forces and stiffnesses
  computeCollisionResponse(R,K,collisionC, true);

  // compute the RHS of the residual:
  TIMER rhsTimer("Forming Initial RHS");
  const REAL invDt = 1.0 / _dt;
  const REAL invDt2 = invDt * invDt;
  _b = (invDt * _M - C) * _velocity + R + _externalForces;

  // collisionC does *not* appear here, since the damping only depends on
  // v^{t+1}, which only depends on \Delta x, which is the variable
  // we are solving for
  //_b = (invDt * _M - C - collisionC) * _velocity + R + _externalForces;
  rhsTimer.stop();

  // assemble system matrix A
  TIMER lhsTimer("Forming Initial LHS");
  // TODO: implement collisions
  //_A = _M * invDt2 - (C + collisionC) * invDt - K; // w/ Rayleigh
  _A = _M * invDt2 - C * invDt - K; // w/ Rayleigh
  lhsTimer.stop();

  // in [TJM15], this is c = b - Az (page 8, top of column 2)
  TIMER projectionTimer("PPCG projection");
  // TODO: add plane constraints
  VECTOR c = _b - _A * z;
  //VECTOR c = _b;

  // just for clarity, go ahead and form the RHS and LHS explicitly
  //
  // Since _S is sparse, this multipy could be accelerated significantly, 
  // but leaving it as it is for now
  VECTOR rhs = _S * c;
  SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
  projectionTimer.stop();

#if 1
  TIMER pcgTimer("PCG Solve");
  _cgSolver.compute(LHS);
  VECTOR y = _cgSolver.solve(rhs);
  pcgTimer.stop();

  if (verbose)
    printf("  PCG iters: %3i err: %6.4e \n", (int)_cgSolver.iterations(), (float)_cgSolver.error());
#else
  TIMER choleskyTimer("Cholesky Solve");
  Eigen::SimplicialLDLT<SPARSE_MATRIX> solver;
  solver.compute(LHS);
  VECTOR y = solver.solve(rhs);
  choleskyTimer.stop();
#endif
  
  // aliasing _solution to \Delta x just to make clear what we're doing here
  VECTOR& xDelta = _solution;
  // TODO: add plane constraints
  xDelta = y + z;
  //xDelta = y;

  // update positions
  _position += xDelta;
  /*
  cout << " x delta: " << xDelta.norm() << endl;
  cout << " b        " << _b.norm() << endl;
  cout << " R        " << R.norm() << endl;
  cout << " ext      " << _externalForces.norm() << endl;
  cout << " Cv:      " << (-C * _velocity).norm() << endl;
  cout << " Ma:      " << (invDt * _M * _velocity).norm() << endl;
  cout << " v:       " << _velocity.norm() << endl;
  cout << " C:       " << C.norm() << endl;
  */

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
    buildConstraintMatrix();
    updateConstraintTargets();
  }
  // update the targets, but the constraint matrix should not have changed.
  else
  {
    updateSurfaceConstraints();
    updateConstraintTargets();
  }
  
  // update node positions
  _triangleMesh.setDisplacement(_position);

  // update velocity
  _velocity = invDt * (_position - _positionOld);

	// update acceleration
  //_acceleration = invDt * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the damping matrix based on the rest pose stiffness
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildRayleighDampingMatrix()
{
  // back up current state
  _temp = _triangleMesh.getDisplacement();

  // set to zero displacement
  VECTOR zero(_DOFs);
  zero.setZero();
  _triangleMesh.setDisplacement(zero);

  // get stiffness matrix at that state
  _triangleMesh.computeFs();
  _triangleMesh.computeSVDs();
  SPARSE_MATRIX K = _triangleMesh.computeStretchingClampedHessian(_stretchingEnergy);

  // restore old state
  _triangleMesh.setDisplacement(_temp);
  
  // build out the Rayleigh damping
  SPARSE_MATRIX C = _rayleighAlpha * _M;
  C += _rayleighBeta * K;

  return C;
}

#if 0
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void printEntry(const VECTOR& v, const int i, const string& varname)
{
  VECTOR3 v3;
  v3[0] = v[3 * i];
  v3[1] = v[3 * i + 1];
  v3[2] = v[3 * i + 2];
  cout << varname.c_str() << ": " << v3.transpose() << endl;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////
// find the surface constraints that are separating
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::findSeparatingSurfaceConstraints(const VECTOR& unfiltered)
{
  bool changed = false;

  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // is the vertex still inside the object? In that case, keep the constraint.
    const KINEMATIC_SHAPE* shape = constraint.shape;
    const int vertexID = constraint.vertexID;
    const REAL signedDistance = shape->signedDistance(_triangleMesh.vertices()[vertexID]);

    bool debug = constraint.vertexID == 262;
    //bool debug = false;
    if (debug)
    {
      cout << " Constraint for vertex: " << constraint.vertexID << endl;
      cout << "   Signed distance:     " << signedDistance << endl;
    }

    /*
    // if the distance is outside and large, move on
    //if (signedDistance > 1e-6) 
    if (signedDistance > 1e-5) 
    //if (signedDistance > 1e-4) 
    {
      constraint.isSeparating = true;
      changed = true;
      if (debug) cout << " CONSTRAINT IS OUTSIDE " << endl;
      continue;
    }
    */

    // what direction is the solution pointing in?
    const int vectorID = 3 * vertexID;
    VECTOR3 xDirection;
    xDirection[0] = unfiltered[vectorID];
    xDirection[1] = unfiltered[vectorID + 1];
    xDirection[2] = unfiltered[vectorID + 2];

    // make the test material agnostic and only look at the direction; then if it's a big force, 
    // the testing threshold won't get messed up later 
    if (xDirection.norm() > 1.0)
      xDirection.normalize();

    // what direction is the kinematic object's surface normal pointing in?
    VECTOR3 normal = shape->localNormalToWorld(constraint.localNormal);

    // what is the magnitude in the separation direction?
    const REAL separationMagnitude = xDirection.dot(normal);

    if (debug) cout << "  separation magnitude: " << separationMagnitude << endl;

    // if the distance is outside and large, move on
    //if (signedDistance > 1e-4 || (signedDistance > 1e-6 && separationMagnitude > 1e-6)) 
    if (signedDistance > 1e-6 && separationMagnitude > 1e-6)
    //if (signedDistance > 1e-5) 
    //if (signedDistance > 1e-4) 
    {
      constraint.isSeparating = true;
      changed = true;
      if (debug) cout << " CONSTRAINT IS OUTSIDE " << endl;
      continue;
    }

    /*
    // is the velocity pulling away?
    VECTOR3 localVelocity = velocity(vertexID);
    const REAL velocitySeparation = localVelocity.dot(normal);
    */

    // are they the same? in that case, they're separating
    //
    // using a small epsilon threshold, especially for velocitySeparation
    // because otherwise there is jittering during resting contact
    if (separationMagnitude > 1e-6)   // velocity condition doesn't seem to ever get triggered
    //if (separationMagnitude > 1e-6 || velocitySeparation > 1e-6)
    //if (separationMagnitude > 1e-5 || velocitySeparation > 1e-5)
    //if (separationMagnitude > 0.0 || velocitySeparation > 0.0)
    //if (separationMagnitude > 1e-4 || velocitySeparation > 1e-4)  // jitter goes away for SNH lambda = 1000, but gets sticky
    //if (separationMagnitude > 9e-5 || velocitySeparation > 9e-5)
    {
      //cout << " Separation: " << separationMagnitude << " velocity: " << velocitySeparation << endl;
      constraint.isSeparating = true;
      changed = true;

      if (debug)
        cout << " SEPARATING" << endl;
    }
  }

  return changed;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// add a gravity body force to the simulation
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::addGravity(const VECTOR3& bodyForce)
{
  const vector<REAL>& oneRingAreas = _triangleMesh.restOneRingAreas();

  for (int x = 0; x < _DOFs / 3; x++)
  {
    const VECTOR3 scaledForce = oneRingAreas[x] * bodyForce;
    _externalForces[3 * x]     += scaledForce[0];
    _externalForces[3 * x + 1] += scaledForce[1];
    _externalForces[3 * x + 2] += scaledForce[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// constrain surface nodes inside a kinematic body to move along with that body
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::attachKinematicSurfaceConstraints(const KINEMATIC_SHAPE* shape)
{
  // get all the nodes inside the shape
  const vector<VECTOR3>& vertices = _triangleMesh.vertices();
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    VECTOR3 v = vertices[x];

    // if it's not inside, move on
    if (!shape->inside(v)) continue;

    // if it's inside, get its local coordinates
    VECTOR3 local = shape->worldVertexToLocal(v);

    // record everything the solver will need later
    KINEMATIC_CONSTRAINT constraint;
    constraint.shape = shape;
    constraint.vertexID = x;
    constraint.localPosition = local;

    // rememeber the constraint for later
    _kinematicConstraints.push_back(constraint);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// release surface nodes inside a kinematic body to move freely
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::releaseKinematicSurfaceConstraints(const KINEMATIC_SHAPE* shape)
{
  // look through all of the constraints and delete the ones associated with the shape
  // inefficient: makes a new list and replaces old _kinematicConstraints.
  vector<KINEMATIC_CONSTRAINT> prunedConstraints;
  for (unsigned int i = 0; i < _kinematicConstraints.size(); i++)
  {
    KINEMATIC_CONSTRAINT c = _kinematicConstraints[i];
    const KINEMATIC_SHAPE* s = c.shape;
    if (shape == s) continue;
    prunedConstraints.push_back(c);
  } 
  _kinematicConstraints = prunedConstraints;

}

#if 0
///////////////////////////////////////////////////////////////////////////////////////////////////////
// constrain all nodes inside a kinematic body to move along with that body
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::attachKinematicConstraints(const KINEMATIC_SHAPE* shape)
{
  // get all the nodes inside the shape
  const vector<VECTOR3>& vertices = _tetMesh.vertices();
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    VECTOR3 v = vertices[x];

    // if it's not inside, move on
    if (!shape->inside(v)) continue;

    // if it's inside, get its local coordinates
    VECTOR3 local = shape->worldVertexToLocal(v);

    // record everything the solver will need later
    KINEMATIC_CONSTRAINT constraint;
    constraint.shape = shape;
    constraint.vertexID = x;
    constraint.localPosition = local;

    // rememeber the constraint for later
    _kinematicConstraints.push_back(constraint);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// which nodes are the constrained ones?
///////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> TIMESTEPPER::constrainedNodes() const
{
  // find the (unique) constrained nodes
  map<int, bool> isConstrained;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const int vertexID = _kinematicConstraints[x].vertexID;
    isConstrained[vertexID] = true;
  }

  // tape out the unique IDs
  vector<int> nodes;
  for (auto iter = isConstrained.begin(); iter != isConstrained.end(); iter++)
    nodes.push_back(iter->first);

  return nodes;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////
// add kinematic collision object to system
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::addKinematicCollisionObject(const KINEMATIC_SHAPE* shape)
{
  // make sure we didn't already add it
  for (unsigned int x = 0; x < _collisionObjects.size(); x++)
    if (shape == _collisionObjects[x])
      return;

  _collisionObjects.push_back(shape);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// find all the surface vertices that are in collision and create constraints
//
// this one is slightly different in QUASISTATIC, i.e. that one can't take
// velocity into account as a separation condition, so this function has
// been made virtual
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::findNewSurfaceConstraints(const bool verbose)
{
  const vector<VECTOR3> vertices = _triangleMesh.vertices();

  if (verbose)
    cout << " Currently tracking " << _planeConstraints.size() << " constraints " << endl;

  // build any new constraints
  int newConstraints = 0;
  for (unsigned int y = 0; y < _collisionObjects.size(); y++)
  {
    const KINEMATIC_SHAPE* shape = _collisionObjects[y];
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // get the vertex
      int vertexID = x;

      //bool debug = (vertexID == 0);
      bool debug = false;

      // if it's already in collision, skip it
      if (_inCollision[vertexID]) 
      {
        if (debug) cout << " vertex is already collision, moving on" << endl;
        continue;
      }

      // see if it's inside the shape
      const VECTOR3& vertex = vertices[vertexID];
      if (!shape->inside(vertex)) 
      {
        if (debug) cout << " vertex is not inside the shape, moving on" << endl;
        continue;
      }

      VECTOR3 closestPoint;
      VECTOR3 closestNormal;
      shape->getClosestPoint(vertex, closestPoint, closestNormal);
     
      // if the velocity is pulling away from the surface, don't constrain it
      VECTOR3 vertexVelocity = velocity(vertexID);
      VECTOR3 normal = shape->localNormalToWorld(closestNormal);
      const REAL velocitySeparation = vertexVelocity.dot(normal);

      if (debug)
      {
        cout << " velocity:   " << vertexVelocity.transpose() << endl;
        cout << " normal:     " << normal.transpose() << endl;
        cout << " separation: " << velocitySeparation << endl;
      }

      // comparing directly against zero here. Trying for a small
      // epsilon just induces sticking.
      //
      // Without this, objects will always stick to a surface after initially
      // sliding
      //
      // BDF-2 sticks unless -FLT_EPSILON is used, but other integrators seem okay
      if (velocitySeparation >= 0)
      //if (velocitySeparation >= -1e-9)
      //if (velocitySeparation >= -FLT_EPSILON)
        continue;

      // store the constraint
      PLANE_CONSTRAINT constraint;
      constraint.shape = shape;
      constraint.vertexID = vertexID;
      constraint.localClosestPoint = closestPoint;
      constraint.localNormal = closestNormal;
      constraint.isSeparating = false;
      addPlaneConstraint(constraint);

      _inCollision[vertexID] = true;
      newConstraints++;
    }
  }
  if (verbose)
    cout << " Found " << newConstraints << " new constraints " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the closest point positions on the surface constraints
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateSurfaceConstraints()
{
  const vector<VECTOR3> vertices = _triangleMesh.vertices();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const KINEMATIC_SHAPE& shape = *_planeConstraints[x].shape;

    // get the new vertex position
    const int vertexID = _planeConstraints[x].vertexID;
    const VECTOR3& vertex = vertices[vertexID];

    // recompute closes point
    VECTOR3 closestPointLocal, normalLocal;
    shape.getClosestPoint(vertex, closestPointLocal, normalLocal);

    // store result
    _planeConstraints[x].localClosestPoint = closestPointLocal;
    _planeConstraints[x].localNormal = normalLocal;
  }
  // we're not checking whether it's still inside or separating here.
  // That will be handled by findSeparatingSurfaceConstraints
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// find all the constraints tagged for deletion and delete them
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::deleteSurfaceConstraints(const bool verbose)
{
  int totalDeleted = 0;

  // If any constraints were tagged for deletion last time, delete them now.
  //
  // That's right, I'm just building a whole new vector instead of deleting nodes 
  // from a linked list. If it's all too ugly for you, look away.
  // Like I said in the README.md, this library is not optimized yet.
  vector<PLANE_CONSTRAINT> constraints;
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    // if it's not separating, just store it
    if (!_planeConstraints[x].isSeparating)
      constraints.push_back(_planeConstraints[x]);
    // if we're deleting this, make sure this surface vertex isn't 
    // tagged as in collision anymore
    else
    {
      _inCollision[_planeConstraints[x].vertexID] = false;
      totalDeleted++;
    }
  }

  if (verbose)
    cout << " Total deleted: " << totalDeleted << endl;

  _planeConstraints = constraints;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// velocity at a specific vertex
///////////////////////////////////////////////////////////////////////////////////////////////////////
const VECTOR3 TIMESTEPPER::velocity(unsigned int index) const
{
  assert(index >= 0);
  assert(index < _velocity.size());
  VECTOR3 vertexVelocity;
  vertexVelocity[0] = _velocity[3 * index];
  vertexVelocity[1] = _velocity[3 * index + 1];
  vertexVelocity[2] = _velocity[3 * index + 2];

  return vertexVelocity;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// reset the Rayleigh damping constants
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::setRayeligh(const REAL alpha, const REAL beta)
{
  _rayleighAlpha = alpha;
  _rayleighBeta = beta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// do the collision detection, in anticipation of collision response
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeCollisionDetection()
{
  TIMER functionTimer(__FUNCTION__);

  // if the tet mesh has an AABB accelerator, refit it
  // TET_MESH_FASTER* fast = dynamic_cast<TET_MESH_FASTER*>(&_tetMesh);
  // if (fast != NULL)
  //   fast->refitAABB();

  // do the collision processing
  const REAL invDt = 1.0 / _dt;
  if (_vertexFaceSelfCollisionsOn)
  {
    // vertex-face collision detection
    _triangleMesh.computeVertexFaceCollisions();

    // build out the vertex-face "collision tets"
    // TODO: this need to get cut down
    _triangleMesh.buildVertexFaceCollisionTets(_velocity);
  }
  if (_edgeEdgeSelfCollisionsOn)
  {
    // edge-edge collision detection
    _triangleMesh.computeEdgeEdgeCollisions();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute collision forces, add them to the forces and stiffness matrix
// R = forces, K = stiffness, C = damping
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeCollisionResponse(VECTOR& R, SPARSE_MATRIX& K, SPARSE_MATRIX& collisionC, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  // build the collision forces and Hessians
  const int rank = R.size();
  VECTOR collisionForces(rank);
  SPARSE_MATRIX collisionK(rank, rank);
  //SPARSE_MATRIX collisionC(rank, rank);

  collisionForces.setZero();
  collisionK.setZero();
  collisionC.setZero();

  const REAL dampingBeta = _collisionDampingBeta;

  // vertex-face case
  VECTOR forcesVF;
  SPARSE_MATRIX hessianVF;
  if (_vertexFaceSelfCollisionsOn)
  {
    // get vertex-face collision forces and gradient
    forcesVF = _triangleMesh.computeVertexFaceCollisionForces(_vfEnergy);
    hessianVF = _triangleMesh.computeVertexFaceCollisionClampedHessian(_vfEnergy);
    
    collisionForces += forcesVF;
    collisionK += hessianVF;
    collisionC += dampingBeta * hessianVF;
  }

  // edge-edge case
  VECTOR forcesEE;
  SPARSE_MATRIX hessianEE;
  if (_edgeEdgeSelfCollisionsOn)
  {

    forcesEE = _triangleMesh.computeEdgeEdgeCollisionForces(_eeEnergy);
    hessianEE = _triangleMesh.computeEdgeEdgeCollisionClampedHessian(_eeEnergy);
    
    collisionForces += forcesEE;
    collisionK += hessianEE;
    collisionC += dampingBeta * hessianEE;
  }

#if VERY_VERBOSE
  if (verbose)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    if (_vertexFaceSelfCollisionsOn)
    {
      cout << " VF f: " << forcesVF.norm() << endl;
      cout << " VF K: " << hessianVF.norm() << endl;
      cout << " VF C: " << (dampingBeta * hessianVF).norm() << endl;
    }

    if (_edgeEdgeSelfCollisionsOn)
    {
      cout << " EE f: " << forcesEE.norm() << endl;
      cout << " EE K: " << hessianEE.norm() << endl;
      cout << " EE C: " << (dampingBeta * hessianEE).norm() << endl;
    }

    cout << " collision forces: " << collisionForces.norm() << endl;
    cout << " collision K: " << collisionK.norm() << endl;
    cout << " collision C: " << collisionC.norm() << endl;
  }
#endif

  // add self-collisions to both LHS and RHS
  if (_vertexFaceSelfCollisionsOn || _edgeEdgeSelfCollisionsOn)
  {
    R += collisionForces;
    K += collisionK;
  }
}

} // HOBAK
} // SHELL
