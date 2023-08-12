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
#ifndef STRAND_VOLUME_HYBRID_TIMESTEPPER_H
#define STRAND_VOLUME_HYBRID_TIMESTEPPER_H

#include "Geometry/STRAND_MESH.h"
#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Strand/STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Volume/EDGE_COLLISION.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Damping/Volume/DAMPING.h"
#include "util/BLOCK_DIAGONAL_MATRIX3.h"
#include "Collision/C_PLANES_EE.h"

#include <iostream>

namespace HOBAK {
namespace STRAND_VOLUME {

////////////////////////////////////////////////////////////////////////////////////////////////////
// This implements Baraff-Witkin-style constraints from
//  "Large Steps in Cloth Simulation", SIGGRAPH 1998
// by building the system described in 
//  "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class TIMESTEPPER
{
public:
  TIMESTEPPER(STRAND_MESH& strandMesh, STRAND::STRETCHING& stretching, TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);
  TIMESTEPPER(STRAND_MESH& strandMesh, STRAND::STRETCHING& stretching, TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, VOLUME::DAMPING& damping);
  virtual ~TIMESTEPPER();

  VECTOR& externalForces()             { return _externalForces; };
  const VECTOR& externalForces() const { return _externalForces; };
  const vector<PLANE_CONSTRAINT>& planeConstraints() const { return _planeConstraints; };
  const vector<KINEMATIC_CONSTRAINT>& kinematicConstraints() const { return _kinematicConstraints; };
  const STRAND::STRETCHING& strand_material() const            { return _stretchingEnergy; };
  const VOLUME::HYPERELASTIC& volume_material() const             { return _hyperelastic; };
  const string strand_materialName() const                     { return _stretchingEnergy.name(); };
  const string volume_materialName() const                        { return _hyperelastic.name(); };
  const string& name() const                            { return _name; };
  const REAL dt() const                                 { return _dt; };
  const REAL rayleighAlphaStrand() const                      { return _rayleighAlphaStrand; };
  const REAL rayleighBetaStrand() const                       { return _rayleighBetaStrand; };
  const REAL rayleighAlphaVolume() const                      { return _rayleighAlphaVolume; };
  const REAL rayleighBetaVolume() const                       { return _rayleighBetaVolume; };

  // todo: add seperate position for strands and volumes 
  const TET_MESH& tetMesh() const       { return _tetMesh; };
  const STRAND_MESH& strandMesh() const       { return _strandMesh; };
  const VECTOR position() const         { return _position; };
  const VECTOR positionOld() const      { return _positionOld; };
  const VECTOR velocity() const         { return _velocity; };
  const VECTOR velocityOld() const { return _velocityOld; };
  VECTOR& velocityOld()            { return _velocityOld; };
  VECTOR& position()                    { return _position; };
  VECTOR& positionOld()                 { return _positionOld; };
  VECTOR& velocity()                    { return _velocity; };
  // Newmark needs to recompute things if this is set differently
  //REAL& dt()                            { return _dt; };
  //
  const bool& vertexFaceSelfCollisionsOn() const { return _vertexFaceSelfCollisionsOn; };
  const bool& edgeEdgeSelfCollisionsOn() const   { return _edgeEdgeSelfCollisionsOn; };
  const REAL& collisionStiffnessStrand() const         { return _collisionStiffnessStrand; };
  const REAL& collisionStiffnessStrandVolume() const         { return _collisionStiffnessStrandVolume; };
  const REAL& collisionDampingBetaStrand() const       { return _collisionDampingBetaStrand; };
  const REAL& collisionDampingBetaStrandVolume() const       { return _collisionDampingBetaStrandVolume; };
  bool& vertexFaceSelfCollisionsOn()    { return _vertexFaceSelfCollisionsOn; };
  bool& edgeEdgeSelfCollisionsOn()      { return _edgeEdgeSelfCollisionsOn; };
  REAL& collisionStiffnessStrand()            { return _collisionStiffnessStrand; };
  REAL& collisionStiffnessStrandVolume()            { return _collisionStiffnessStrandVolume; };
  REAL& collisionDampingBetaStrand()          { return _collisionDampingBetaStrand; };
  REAL& collisionDampingBetaStrandVolume()          { return _collisionDampingBetaStrandVolume; };
  const REAL& collisionStiffnessVolume() const         { return _collisionStiffnessVolume; };
  const REAL& collisionDampingBetaVolume() const       { return _collisionDampingBetaVolume; };
  REAL& collisionStiffnessVolume()            { return _collisionStiffnessVolume; };
  REAL& collisionDampingBetaVolume()          { return _collisionDampingBetaVolume; };
  virtual void setDt(const REAL dt)     { _dt = dt; };
  void setRayelighStrand(const REAL alpha, const REAL beta);
  void setRayelighVolume(const REAL alpha, const REAL beta);

  // feature toggles
  bool& collisionsEnabled()       { return _collisionsEnabled; };
  bool& pcgEnabled()              { return _pcgEnabled; };
  bool& hessianClampingEnabled()  { return _hessianClampingEnabled; };
  unsigned int& maxNewtonIterations() { return _maxNewtonIterations; };
  const unsigned int& maxNewtonIterations() const { return _maxNewtonIterations; };
  bool& disablePreconditioner()  { return _disablePreconditioner; };

  // velocity at a specific vertex
  const VECTOR3 velocity(unsigned int index) const;

  // take a timestep
  virtual bool solveDynamics(const bool verbose);
  bool solveEnergyDamped(const bool verbose);
  // solves with self collisions, and collision forces are Rayleigh damped
  bool solveRayleighDamped(const bool verbose);

  // this isn't quite quasistatic, but will do for a unit test
  bool solveQuasistatic(const bool verbose);

  // take a timestep, with multiple Newton iterations
  bool solveNewton(const bool verbose);

  // add a gravity body force to the simulation
  void addGravity(const VECTOR3& bodyForce);

  // add a plane constraint
  void addPlaneConstraint(const PLANE_CONSTRAINT& constraint) { _planeConstraints.push_back(constraint); };
  void clearPlaneConstraints()                                { _planeConstraints.clear(); };
  int totalPlaneConstraints()                                 { return _planeConstraints.size(); };

  // constrain surface nodes inside a kinematic body to move along with that body
  void attachKinematicSurfaceConstraints(const KINEMATIC_SHAPE* shape, const bool verbose = false);
  void removeKinematicSurfaceConstraints();
  void popBackKinematicSurfaceConstraints();

  // constrain all nodes inside a kinematic body to move along with that body
  void attachKinematicConstraints(const KINEMATIC_SHAPE* shape);

  // which nodes are the constrained ones?
  vector<int> constrainedNodes() const;

  // add kinematic collision object to system
  void addKinematicCollisionObject(const KINEMATIC_SHAPE* shape);

  // make all objects lighter or heavier
  void scaleMass(const REAL& scalar)  { _M *= scalar; };

  // apply a rigid rotation as a warm start to the solve
  void applyRigidRotation(const MATRIX3& R, const VECTOR3& translation);

  // see if a specific vertex is kinematically constrained
  bool isKinematicallyConstrained(const int vertexID);

protected:
  // shared initialization across constructors
  void initialize();

  // build the constraint matrix to incorporate Baraff-Witkin-style constraints, 
  // but using the [TJM15] projection matrix
  virtual void buildConstraintMatrix();
  void buildConstraintMatrixEdgeEnd();
  void buildConstraintMatrixInterleaved();
  virtual void buildConstraintMatrixFaster();
  void buildConstraintMatrixFasterEdgeEnd();
  void buildConstraintMatrixFasterInterleaved();
  void buildBlockConstraintMatrix();
  SPARSE_MATRIX& filteredSystem();
  SPARSE_MATRIX buildPlaneConstraintsNoIdentity();
  SPARSE_MATRIX buildPlaneConstraintsOnly();

  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  //
  // the target formulations are different in dynamics vs. quasistatics, and
  // position vs. velocity updates, so it is pure virtual here
  //virtual void updateConstraintTargets() = 0;
  virtual void updateConstraintTargets();
  void updateConstraintTargetsEdgeEnd();
  void updateConstraintTargetsInterleaved();

  // filter positions to incorporate Baraff-Witkin-style constraints,
  // this one is slightly different in QUASISTATIC, so this functions has been
  // made virtual
  virtual void applyKinematicConstraints();
  void applyKinematicConstraintsEdgeEnd();
  void applyKinematicConstraintsInterleaved();

  // find all the surface vertices that are in collision and create constraints
  // this one is slightly different in QUASISTATIC, i.e. that one can't take
  // velocity into account as a separation condition, so this function has
  // been made virtual
  virtual void findNewSurfaceConstraints(const bool verbose = false);

  // update the closest point positions on surface constraints
  void updateSurfaceConstraints();

  // find all the constraints tagged for deletion and delete them
  void deleteSurfaceConstraints(const bool verbose = false);

  // find the surface constraints that are separating
  bool findSeparatingSurfaceConstraints(const VECTOR& unfiltered);

  // build the mass matrix based on the one-ring volumes
  virtual SPARSE_MATRIX buildMassMatrix();
  SPARSE_MATRIX buildMassMatrixEdgeEnd();
  SPARSE_MATRIX buildMassMatrixInterleaved();
  void buildMassMatrixVolume(vector<Eigen::Triplet<REAL>>& triplets);
  
  // build the damping matrix based on the rest pose stiffness
  virtual SPARSE_MATRIX buildRayleighDampingMatrix();

  // do the collision detection, in anticipation of collision response
  void computeCollisionDetectionVolume();
  void computeCollisionDetectionStrandVolume();

  // compute collision forces, add them to the forces and stiffness matrix
  // R = forces, K = stiffness, C = damping
  void computeCollisionResponse(VECTOR& R, SPARSE_MATRIX& K, SPARSE_MATRIX& C, const bool verbose = false);
  VECTOR computeEdgeEdgeCollisionForcesStrandVolume() const;
  SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessianStrandVolume() const;
  VECTOR computeVertexFaceCollisionForcesStrandVolume() const;
  SPARSE_MATRIX computeVertexFaceCollisionClampedHessianStrandVolume() const;
  // once attachKinematicSurfaceConstraints is called, we can also see which edges
  // should be constrained
  void computeConstrainedEdges(const bool verbose = false);

  // indirect function so we can try out lots of different solvers
  virtual VECTOR solveSystem(const SPARSE_MATRIX& LHS, const VECTOR& rhs);

  // build out the PPCG projection
  SPARSE_MATRIX constraintProjectFaster();

  // indexing
  const VECTOR3& getMeshVertex(const int vertexID) const;
  const VECTOR3& getRestMeshVertex(const int vertexID) const;
  const int getPositionIndexEdgeEnd(const int vertexID) const;
  const int getPositionIndexInterleaved(const int vertexID) const;
  void concatSparseMatrix(const SPARSE_MATRIX& mStrand, const SPARSE_MATRIX& mVolume, SPARSE_MATRIX& mat);

  REAL _residual;
  int _seenPCGIterations;

  STRAND_MESH& _strandMesh;
  STRAND::STRETCHING& _stretchingEnergy;
  //STRAND::BENDING& _bendingEnergy;
  //VOLUME::DAMPING* _damping;

  TET_MESH& _tetMesh;
  VOLUME::HYPERELASTIC& _hyperelastic;
  VOLUME::DAMPING* _damping;

  int _DOFs;
  int _DOFsStrand;
  int _DOFsVolume;
  VECTOR _forces;
  VECTOR _externalForces;
  int _totalVertices;

  // RHS of the solve
  VECTOR _b;

  // result of the solve
  VECTOR _solution;

  // everybody needs a scratchpad sometimes
  VECTOR _temp;

  // constraint matrix
  SPARSE_MATRIX _S;
  SPARSE_MATRIX _IminusS;
  SPARSE_MATRIX _I;
  SPARSE_MATRIX _N;
  // build block versions of the constraint matrices
  BLOCK_DIAGONAL_MATRIX3 _blockS;
  BLOCK_DIAGONAL_MATRIX3 _blockIminusS;

  // constraint targets
  VECTOR _constraintTargets;

  // is the vertex already experiencing a kinematic collision?
  vector<bool> _inCollision;

  // constraints to have vertices move with a kinematic body
  vector<KINEMATIC_CONSTRAINT> _kinematicConstraints;

  // constraints to have vertex slide along a plane
  vector<PLANE_CONSTRAINT> _planeConstraints;

  // kinematic collision objects
  vector<const KINEMATIC_SHAPE*> _collisionObjects;

  // list of constrained edges
  vector<int> _constrainedEdges;

  // variables to solve for
  VECTOR _position;
  VECTOR _velocity;
  VECTOR _velocityDelta;

  // in case the user wants to rewind to the previous positions 
  VECTOR _positionOld;

  // timestep
  REAL _dt;
  REAL _rayleighAlphaStrand;
  REAL _rayleighBetaStrand;
  REAL _rayleighAlphaVolume;
  REAL _rayleighBetaVolume;
  
  // solver vars
  SPARSE_MATRIX _A;
  SPARSE_MATRIX _M;
  SPARSE_MATRIX _MStrand;
  SPARSE_MATRIX _MVolume;
  
  // global Hessian matrix
  SPARSE_MATRIX _H;
  Eigen::ConjugateGradient<SPARSE_MATRIX, Eigen::Lower|Eigen::Upper> _cgSolver;

  // what's this timestepper called?
  string _name;

  // are self-collisions activated?
  bool _vertexFaceSelfCollisionsOn;
  bool _edgeEdgeSelfCollisionsOn;

  // collision spring and damping constants
  REAL _collisionStiffnessStrand;
  REAL _collisionDampingBetaStrand;
  REAL _collisionStiffnessVolume;
  REAL _collisionDampingBetaVolume;
  REAL _collisionStiffnessStrandVolume;
  REAL _collisionDampingBetaStrandVolume;

  // collisions between strands and volumes
  // first ID into strand vertex 
  // seconde ID into _tetmesh._surfaceTranigles
  vector<pair<int, int> > _vertexFaceCollisionsStrandVolume;
  // list of "collision tets" formed by vertex-face pairs
  // first vertex on strand. the other 3 on volume
  vector<VECTOR4I> _vertexFaceCollisionTetsStrandVolume;
  // list of "collision tets" formed by edge-edge pairs
  //vector<VECTOR4I> _edgeEdgeCollisionsTets;
  // scaling term for vertex-face collision forces
  vector<REAL> _vertexFaceCollisionAreasStrandVolume;

  // list of edge-edge collision indices
  // first indexes into _surfaceEdges
  // second indexes into _strandMesh.edgeIndices()
  vector<pair<int, int> > _edgeEdgeCollisionsStrandVolume;
  // interpolation coordinates for edge-edge collisions
  vector<pair<VECTOR2, VECTOR2> > _edgeEdgeCoordinatesStrandVolume;
  // are the edge-edge collisions still separate, or is there already a face-edge intersection?
  vector<bool> _edgeEdgeIntersectionsStrandVolume;
  // scaling term for edge-edge collision forces
  vector<REAL> _edgeEdgeCollisionAreasStrandVolume;

  // Newton iteration variables
  REAL _residualTolerance;
  unsigned int _maxNewtonIterations;

  // BDF-1 variables
  REAL _time;
  int _currentTimestep;
  VECTOR _velocityOld;
  VECTOR _acceleration;

  // toggle various features
  bool _collisionsEnabled;
  bool _pcgEnabled;
  bool _hessianClampingEnabled;
  bool _disablePreconditioner;

  // preallocated sparse matrices for fast constraint projection 
  SPARSE_MATRIX _BA, _C, _P, _PnoI, _PAP;

  // preallocated sparse matrices for missing term in fast constraint projection 
  SPARSE_MATRIX _missing, _blockDiag;

  // which vertex-face collision force are we using?
  VOLUME::VERTEX_FACE_COLLISION* _vertexFaceEnergyStrandVolume;
  // which edge-edge collision force are we using?
  VOLUME::EDGE_COLLISION* _edgeEdgeEnergyStrandVolume;
  REAL _collisionEps;
};

} // HOBAK
} // TIMESTEPPER

#endif
