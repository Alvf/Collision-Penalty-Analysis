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
#include "Hyperelastic/Volume/VERTEX_FACE_COLLISION.h"
#include "Hyperelastic/Volume/VERTEX_FACE_SQRT_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_SQRT_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_HYBRID_COLLISION.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "Damping/Volume/ARAP_DAMPING.h"
#include "Geometry/TET_MESH_FASTER.h"
#include "Geometry/TET_STRAND_MESH.h"
#include "TIMER.h"
#include <float.h>
#include <iostream>
#include "PCG.h"
#include "DIAGONAL.h"
#include "STRAND_DIAGONAL.h"

#include <Eigen/SparseLU>

#define USING_AMGCL 0

#if USING_AMGCL
/* Example is from:
https://gist.github.com/ddemidov/ac9a060507edea7586d0516b174ccd13#file-eigen-cpp
https://gist.githubusercontent.com/ddemidov/ac9a060507edea7586d0516b174ccd13/raw/e572d888a3b1138d9fa5d68fea6287f245f10168/eigen.cpp
 */
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/backend/eigen.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/preconditioner/dummy.hpp>
#endif

#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

static Eigen::IOFormat octave(-1, 0, ", ", ";\n", "", "", "[", "]");

using namespace std;

namespace HOBAK {
namespace STRAND_VOLUME {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
TIMESTEPPER::TIMESTEPPER(STRAND_MESH& strandMesh, STRAND::STRETCHING& stretching, 
                         TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  _strandMesh(strandMesh), _stretchingEnergy(stretching), 
  _tetMesh(tetMesh), _hyperelastic(hyperelastic), _damping(NULL)
{
  initialize();

  //computeRestTwistAngles();
}

TIMESTEPPER::TIMESTEPPER(STRAND_MESH& strandMesh, STRAND::STRETCHING& stretching, 
                         TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, 
                         VOLUME::DAMPING& damping) :
  _strandMesh(strandMesh), _stretchingEnergy(stretching), 
  _tetMesh(tetMesh), _hyperelastic(hyperelastic), _damping(&damping)
{
  initialize();

  //computeRestTwistAngles();
}

TIMESTEPPER::~TIMESTEPPER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::initialize()
{
  cout << " Initializing strand volume integrator ... "<<endl;

  _residual = FLT_MAX;
  _seenPCGIterations    = -1;
  _currentTimestep = 0;

  _DOFsStrand = _strandMesh.DOFs();
  // cout<<"_DOFsStrand: " << _DOFsStrand << endl;
  _DOFsVolume = _tetMesh.DOFs();
  _DOFs = _DOFsStrand + _DOFsVolume;
  _b.resize(_DOFs);
  _forces.resize(_DOFs);
  _externalForces.resize(_DOFs);
  _constraintTargets.resize(_DOFs);

  _b.setZero();
  _forces.setZero();
  _externalForces.setZero();
  _constraintTargets.setZero();

  _totalVertices = _strandMesh.vertices().size() + _tetMesh.vertices().size();
  _inCollision.resize(_totalVertices);
  for (int x = 0; x < _totalVertices; x++)
    _inCollision[x] = false;

  // strand first volume next
  _position.resize(_DOFs);
  _positionOld.resize(_DOFs);
  _velocity.resize(_DOFs);
  _velocityDelta.resize(_DOFs);
  _temp.resize(_DOFs);
  _solution.resize(_DOFs);
  _acceleration.resize(_DOFs);
  _velocityOld.resize(_DOFs);
  // don't just assume it's zero, extract it from the mesh
  //_position.setZero();
  //_positionOld.setZero();

  VECTOR displacement;
  displacement.resize(_DOFs);
  displacement << _strandMesh.getDisplacement(), _tetMesh.getDisplacement();
  _position = displacement;
  _positionOld = _position;
  _velocity.setZero();
  _velocityDelta.setZero();
  _temp.setZero();
  _solution.setZero();
  _positionOld.setZero();
  _acceleration.setZero();
  _velocityOld.setZero();

  _name = string("STRAND_VOLUME::TIMESTEPPER");

  _vertexFaceSelfCollisionsOn = true;
  _edgeEdgeSelfCollisionsOn   = true;

  _dt = 1.0 / 30.0;
  // _rayleighAlphaStrand = 0.001;
  // _rayleighBetaStrand = 0.001;
  _rayleighAlphaVolume = 0.1;
  _rayleighBetaVolume = 0.1;
  //_rayleighAlpha = 0.001;
  //_rayleighBeta = 0.001;
  _rayleighAlphaStrand = 0.0;
  _rayleighBetaStrand = 0.0;

  // build the mass matrix
  _M = buildMassMatrix();
  _MStrand = _M.block(0, 0, _DOFsStrand, _DOFsStrand);
  _MVolume = _M.block(_DOFsStrand, _DOFsStrand, _DOFsVolume, _DOFsVolume);
  // cout<<"M Strand norm: "<<_MStrand.norm() / _MStrand.rows()/_MStrand.cols()<<endl;
  // cout<<"M Volume norm: "<<_MVolume.norm() / _MVolume.rows()/_MVolume.cols()<<endl;

  
  //_collisionStiffness = 10000.0;
  _collisionStiffnessStrand = 1000;
  //_collisionDampingBeta = 0.0;
  _collisionDampingBetaStrand = 0.001;
  //_collisionDampingBeta = 0.01;
  //_collisionDampingBeta = 0.1;
  _collisionStiffnessVolume = 1.0;
  _collisionDampingBetaVolume = 0.001;
  _collisionStiffnessStrandVolume = 1e5;
  _collisionDampingBetaStrandVolume = 0.001;

  _vertexFaceEnergyStrandVolume = new VOLUME::VERTEX_FACE_SQRT_COLLISION(_collisionStiffnessStrandVolume, _collisionEps); // default
  _edgeEdgeEnergyStrandVolume = new VOLUME::EDGE_SQRT_COLLISION(_collisionStiffnessStrandVolume, _collisionEps); // default

  _residualTolerance = 1e-2;
  //_residualTolerance = 1e-4;
  _maxNewtonIterations = 10;

  // feature defaults
  _collisionsEnabled = true;
  _pcgEnabled = true;
  _hessianClampingEnabled = true;

  _I = SPARSE_MATRIX(_DOFs, _DOFs);
  _I.setIdentity();

  _disablePreconditioner = true;

  cout << "done. " << endl;
}

// indexing
const VECTOR3& TIMESTEPPER::getMeshVertex(const int vertexID) const
{
  if(vertexID >= _strandMesh.totalVertices()){
      const int volumeVID = vertexID - _strandMesh.totalVertices();
      return _tetMesh.vertices()[volumeVID];
  }
  else{
    return _strandMesh.vertices()[vertexID];
  }
}

const VECTOR3& TIMESTEPPER::getRestMeshVertex(const int vertexID) const
{
  VECTOR3 vertex;
  if(vertexID >= _strandMesh.totalVertices()){
      const int volumeVID = vertexID - _strandMesh.totalVertices();
      return _tetMesh.restVertices()[volumeVID];
  }
  else{
    return _strandMesh.restVertices()[vertexID];
  }
}

const int TIMESTEPPER::getPositionIndexEdgeEnd(const int vertexID) const
{
  int index;
  if(vertexID >= _strandMesh.totalVertices()){
      const int volumeVID = vertexID - _strandMesh.totalVertices();
      index = _DOFsStrand + 3 * volumeVID;
  }
  else{
    index = 3 * vertexID;
  }
  return index;
}

const int TIMESTEPPER::getPositionIndexInterleaved(const int vertexID) const
{
  int index;
  const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  if(vertexID >= _strandMesh.totalVertices()){
      const int volumeVID = vertexID - _strandMesh.totalVertices();
      index = _DOFsStrand + 3 * volumeVID;
  }
  else {
    index = globalVertexIndices[vertexID];
  }
  return index;
}

void TIMESTEPPER::concatSparseMatrix(const SPARSE_MATRIX& mStrand, const SPARSE_MATRIX& mVolume, SPARSE_MATRIX& mat) {
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int k = 0; k < mStrand.outerSize(); ++k){
    for(SPARSE_MATRIX::InnerIterator it(mStrand, k); it; ++it) {
      TRIPLET triplet(it.row(), it.col(), it.value());
      triplets.push_back(triplet);
    }
  }
  for (unsigned int k = 0; k < mVolume.outerSize(); ++k){
    for(SPARSE_MATRIX::InnerIterator it(mVolume, k); it; ++it) {
      TRIPLET triplet(it.row() + _DOFsStrand, it.col() + _DOFsStrand, it.value());
      triplets.push_back(triplet);
    }
  }
  mat.setFromTriplets(triplets.begin(), triplets.end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the displacement targets the the Baraff-Witkin-style constraints
// are trying to hit. Assumes that buildConstraintMatrix() has already been called
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateConstraintTargets()
{
  //std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  //cout << " Constraint targets: " << _planeConstraints.size() << endl;
  if (_strandMesh.edgeEnd())
  {
    updateConstraintTargetsEdgeEnd();
    return;
  }
    
  updateConstraintTargetsInterleaved();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateConstraintTargetsEdgeEnd()
{
  // timing is negligible
  //TIMER functionTimer(__FUNCTION__);
  _constraintTargets.setZero();

  vector<bool> isKinematic(_totalVertices);
  for (int x = 0; x < _totalVertices; x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    // should ignore if we've tagged it for deletion
    if (_planeConstraints[x].isSeparating)
      continue;

    // retrieve collision information
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const int vertexID = constraint.vertexID;

    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    // compute the target displacement
    const VECTOR3& vertex = getMeshVertex(vertexID);
    const int index = getPositionIndexEdgeEnd(vertexID);
    const VECTOR3& localClosestPoint = _planeConstraints[x].localClosestPoint;
    const VECTOR3& closestPoint = shape->localVertexToWorld(localClosestPoint);

    const VECTOR3& displacement = closestPoint - vertex;
    for (int i = 0; i < 3; i++)
      _constraintTargets[index + i] = displacement[i];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateConstraintTargetsInterleaved()
{
  //TIMER functionTimer(__FUNCTION__);

  vector<bool> isKinematic(_totalVertices);
  for (int x = 0; x < _totalVertices; x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

  _constraintTargets.setZero();
  // const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    // should ignore if we've tagged it for deletion
    if (_planeConstraints[x].isSeparating)
      continue;

    // retrieve collision information
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const int vertexID = constraint.vertexID;
    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    // compute the target displacement
    const VECTOR3& vertex = getMeshVertex(vertexID);
    // cout<<"vertex: "<< vertex << endl;
    // cout<<"vertexID: "<< vertexID << endl;
    // cout<<"_strandMesh.vertices()[vertexID]: "<< _strandMesh.vertices()[vertexID] << endl;
    const int index = getPositionIndexInterleaved(vertexID);
    const VECTOR3& localClosestPoint = _planeConstraints[x].localClosestPoint;
    const VECTOR3& closestPoint = shape->localVertexToWorld(localClosestPoint);


    const VECTOR3& displacement = closestPoint - vertex;
    // cout<<"displacement: "<< displacement << endl;
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
  if (_strandMesh.edgeEnd())
  {
    applyKinematicConstraintsEdgeEnd();
    return;
  }
  applyKinematicConstraintsInterleaved();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::applyKinematicConstraintsEdgeEnd()
{
  // timing is negligible
  //TIMER functionTimer(__FUNCTION__);
  cout << " Kinematic constraints: " << _kinematicConstraints.size() << endl;
  const vector<VECTOR3>& restVerticesStrand = _strandMesh.restVertices();
  const vector<VECTOR3>& restVerticesVolume = _tetMesh.restVertices();
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // pin the mesh position according to the constraint
    const VECTOR3& localPosition = constraint.localPosition;
    VECTOR3 world = constraint.shape->localVertexToWorld(localPosition);

    const int vertexID = constraint.vertexID;
    const VECTOR3& restVertex = getRestMeshVertex(vertexID);
    const int index = getPositionIndexEdgeEnd(vertexID);
    const VECTOR3 diff = world - restVertex;
    _position[index] = diff[0];
    _position[index + 1] = diff[1];
    _position[index + 2] = diff[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::applyKinematicConstraintsInterleaved()
{
  //TIMER functionTimer(__FUNCTION__);
  const vector<VECTOR3>& restVerticesStrand = _strandMesh.restVertices();
  const vector<VECTOR3>& restVerticesVolume = _tetMesh.restVertices();
  const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // pin the  mesh position according to the constraint
    const VECTOR3& localPosition = constraint.localPosition;
    VECTOR3 world = constraint.shape->localVertexToWorld(localPosition);

    const int vertexID = constraint.vertexID;
    const VECTOR3& restVertex = getRestMeshVertex(vertexID);
    const int index = getPositionIndexInterleaved(vertexID);
    const VECTOR3 diff = world - restVertex;
    _position[index] = diff[0];
    _position[index + 1] = diff[1];
    _position[index + 2] = diff[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the constraint matrix to incorporate Baraff-Witkin-style constraints, but using 
// the "Smoothed aggregation multigrid for cloth simulation" projection matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrix()
{
  if (_strandMesh.edgeEnd())
  {
    buildConstraintMatrixEdgeEnd();
    return;
  }

  buildConstraintMatrixInterleaved();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the constraint matrix to incorporate Baraff-Witkin-style constraints, but using 
// the "Smoothed aggregation multigrid for cloth simulation" projection matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrixEdgeEnd()
{
  TIMER functionTimer(__FUNCTION__);
  SPARSE_MATRIX I(_DOFs, _DOFs);
  I.setIdentity();
  _S = I;

  int totalVertices = _totalVertices;
  vector<bool> isKinematic(totalVertices);
  for (int x = 0; x < totalVertices; x++)
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
    const MATRIX3 Sblock = MATRIX3::Identity() - normal * normal.transpose();
    const int vertexID = constraint.vertexID;

    // if it's kinematic, move on
    if (isKinematic[vertexID]) continue;

    const int index = getPositionIndexEdgeEnd(vertexID);
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
    const int index = getPositionIndexEdgeEnd(constraint.vertexID);
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        _S.coeffRef(index + i, index + j) = 0.0;
  }
 
  // unsigned int vertexEnd = 3 * _strandMesh.totalVertices();
  // for (unsigned int x = 0; x < _constrainedEdges.size(); x++)
  // {
  //   unsigned int index = vertexEnd + _constrainedEdges[x];
  //   _S.coeffRef(index,index) = 0.0;
  // }

  // store the complement
  _IminusS = I - _S;
  _N = _IminusS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the constraint matrix to incorporate Baraff-Witkin-style constraints, but using 
// the "Smoothed aggregation multigrid for cloth simulation" projection matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrixInterleaved()
{
  TIMER functionTimer(__FUNCTION__);
  SPARSE_MATRIX I(_DOFs, _DOFs);
  I.setIdentity();
  _S = I;

  int totalVertices = _totalVertices;
  vector<bool> isKinematic(totalVertices);
  for (int x = 0; x < totalVertices; x++)
    isKinematic[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
    isKinematic[_kinematicConstraints[x].vertexID] = true;

  // build the plane constraints for the LHS
  const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
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

    const int index = getPositionIndexInterleaved(vertexID);

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
    const int index = getPositionIndexInterleaved(constraint.vertexID);
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        _S.coeffRef(index + i, index + j) = 0.0;
  }
 
  // const VECTORI& globalEdgeIndices = _strandMesh.globalEdgeIndices();
  // for (unsigned int x = 0; x < _constrainedEdges.size(); x++)
  // {
  //   unsigned int index = globalEdgeIndices[_constrainedEdges[x]];
  //   _S.coeffRef(index,index) = 0.0;
  // }

  // store the complement
  _IminusS = I - _S;
  _N = _IminusS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the constraint matrix to incorporate Baraff-Witkin-style constraints, but using 
// the "Smoothed aggregation multigrid for cloth simulation" projection matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrixFaster()
{
  if (_strandMesh.edgeEnd())
  {
    buildConstraintMatrixFasterEdgeEnd();
    return;
  }
  buildConstraintMatrixFasterInterleaved();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrixFasterEdgeEnd()
{
  TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsS;
  vector<TRIPLET> tripletsN;

  vector<bool> diagonalSeen(_DOFs);
  for (int x = 0; x < _DOFs; x++)
    diagonalSeen[x] = false;

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
    const int index = getPositionIndexEdgeEnd(vertexID);
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
    const int index = getPositionIndexEdgeEnd(constraint.vertexID);

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

  // zero out constraint rows
  // unsigned int vertexEnd = 3 * _strandMesh.totalVertices();
  // for (unsigned int x = 0; x < _constrainedEdges.size(); x++)
  // {
  //   unsigned int index = vertexEnd + _constrainedEdges[x];
  //   diagonalSeen[index] = true;
  //   TRIPLET tripletS(index, index, 0);
  //   tripletsS.push_back(tripletS);
    
  //   TRIPLET tripletN(index, index, 1);
  //   tripletsN.push_back(tripletN);
  // }

  // if the diagonal was never set, set it to one
  for (int x = 0; x < _DOFs; x++)
  {
    if (diagonalSeen[x]) continue;
    TRIPLET tripletS(x, x, 1);
    tripletsS.push_back(tripletS);
  }

  _S = SPARSE_MATRIX(_DOFs, _DOFs);
  _S.setFromTriplets(tripletsS.begin(), tripletsS.end());
  
  _N = SPARSE_MATRIX(_DOFs, _DOFs);
  _N.setFromTriplets(tripletsN.begin(), tripletsN.end());

  // store the complement
  //_IminusS = _I - _S;
  _IminusS = _N;

  /*
  SPARSE_MATRIX diff = (_N - (_I - _S));
  cout << " I - S Diff: " << diff.norm() << endl;
  diff = _S - (_I - _N);
  cout << " S Diff: " << diff.norm() << endl;
  */
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrixFasterInterleaved()
{
  TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsS;
  vector<TRIPLET> tripletsN;

  vector<bool> diagonalSeen(_DOFs);
  for (int x = 0; x < _DOFs; x++)
    diagonalSeen[x] = false;

  // build the plane constraints for the LHS
  // const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
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
    const int index = getPositionIndexInterleaved(vertexID);
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
    const int index = getPositionIndexInterleaved(constraint.vertexID);

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

  // zero out constraint rows
  // const VECTORI& globalEdgeIndices = _strandMesh.globalEdgeIndices();
  // for (unsigned int x = 0; x < _constrainedEdges.size(); x++)
  // {
  //   unsigned int index = globalEdgeIndices[_constrainedEdges[x]];
  //   diagonalSeen[index] = true;
  //   TRIPLET tripletS(index, index, 0);
  //   tripletsS.push_back(tripletS);
    
  //   TRIPLET tripletN(index, index, 1);
  //   tripletsN.push_back(tripletN);
  // }

  // if the diagonal was never set, set it to one
  for (int x = 0; x < _DOFs; x++)
  {
    if (diagonalSeen[x]) continue;
    TRIPLET tripletS(x, x, 1);
    tripletsS.push_back(tripletS);
  }

  _S = SPARSE_MATRIX(_DOFs, _DOFs);
  _S.setFromTriplets(tripletsS.begin(), tripletsS.end());
  
  _N = SPARSE_MATRIX(_DOFs, _DOFs);
  _N.setFromTriplets(tripletsN.begin(), tripletsN.end());

  // store the complement
  //_IminusS = _I - _S;
  _IminusS = _N;

  /*
  SPARSE_MATRIX diff = (_N - (_I - _S));
  cout << " I - S Diff: " << diff.norm() << endl;
  diff = _S - (_I - _N);
  cout << " S Diff: " << diff.norm() << endl;
  */
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildMassMatrix()
{
  if (_strandMesh.edgeEnd())
    return buildMassMatrixEdgeEnd();
  return buildMassMatrixInterleaved();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the volume mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildMassMatrixVolume(vector<Eigen::Triplet<REAL>>& triplets)
{
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  const vector<REAL>& volumes = _tetMesh.restOneRingVolumes(); 

  // set diagonal the one-ring volumes
  // TODO: take into account density
  for (int x = 0; x < _tetMesh.totalVertices(); x++)
  {
    const REAL entry = volumes[x];
    for (int y = 0; y < 3; y++)
    {
      TRIPLET triplet(_DOFsStrand + 3 * x + y, _DOFsStrand + 3 * x + y, entry);
      triplets.push_back(triplet);
    }
  }
  
  // SPARSE_MATRIX A(_DOFs, _DOFs);
  // A.setFromTriplets(triplets.begin(), triplets.end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the strand mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildMassMatrixEdgeEnd()
{
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;

  // set diagonal the one-ring volumes
  for (int x = 0; x < _strandMesh.totalVertices(); x++)
  {
    const REAL entry = _strandMesh.vertexMass(x);
    for (int y = 0; y < 3; y++)
    {
      TRIPLET triplet(3 * x + y, 3 * x + y, entry);
      triplets.push_back(triplet);
    }
  }

  buildMassMatrixVolume(triplets);
  // edges are implicitly set to zero, so nothing left to do
 
  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the strand mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildMassMatrixInterleaved()
{
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;

  // set diagonal the one-ring volumes
  const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  for (int x = 0; x < _strandMesh.totalVertices(); x++)
  {
    const REAL entry = _strandMesh.vertexMass(x);
    const int index = globalVertexIndices[x];
    for (int y = 0; y < 3; y++)
    {
      TRIPLET triplet(index + y, index + y, entry);
      triplets.push_back(triplet);
    }
  }

  buildMassMatrixVolume(triplets);
 
  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;

  return A;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// this isn't quite quasistatic, but will do for a unit test
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::solveQuasistatic(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " STRAND_VOLUME::TIMESTEPPER RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // build the mass matrix
  _M = buildMassMatrix();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;

  // should need to call once, but then preserved throughout
  //applyKinematicConstraints();

  // store the filtered b for later
  //VECTOR unfiltered;

  // build new constraints and see if we should break any
  //findNewSurfaceConstraints(verbose);
  buildConstraintMatrix();

  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  if (tetStrandMesh == NULL)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " Wrong STRAND_MESH type passed to VOLUME_TIMESTEPPER!!!" << std::endl;
    assert(tetStrandMesh);
  }
  tetStrandMesh->computeFs();
  tetStrandMesh->computeSVDs();

  //_strandMesh.setDisplacement(_position);

  // TODO: add plane constraints so that this starts firing
  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  //updateConstraintTargets();
  VECTOR z = _IminusS * _constraintTargets;

  // get the internal forces
  const VECTOR stretchingForce = tetStrandMesh->computeStretchingForces();
  const VECTOR bendingForce = tetStrandMesh->computeBendingForces();
  const VECTOR twistingForce = tetStrandMesh->computeTwistingForces();

  const VECTOR RStrand = bendingForce + stretchingForce + twistingForce;
  const VECTOR RVolume = _tetMesh.computeInternalForces(_hyperelastic, *_damping);
  VECTOR R;
  R.resize(_DOFs);
  R << RStrand, RVolume;
  //std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  //std::cout << " Internal forces: " << std::endl;
  //std::cout << "mine = " << R.transpose().format(octave) << std::endl;
  //std::cout << "stretchingForceMine = " << stretchingForce.transpose().format(octave) << std::endl;
  //std::cout << "bendingForceMine = " << bendingForce.transpose().format(octave) << std::endl;
  //std::cout << "twistingForceMine = " << twistingForce.transpose().format(octave) << std::endl;

  // get the reduced stiffness matrix
  SPARSE_MATRIX KStrand = tetStrandMesh->computeHyperelasticClampedHessianFast();
  // KStrand += _strandMesh.computeBendingHessian();
  // KStrand += _strandMesh.computeTwistingHessian();
  SPARSE_MATRIX KVolume = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

  SPARSE_MATRIX K(_DOFs, _DOFs);
  concatSparseMatrix(KStrand, KVolume, K);

  //SPARSE_MATRIX K = _strandMesh.computeStretchingClampedHessian();
  //K += _strandMesh.computeBendingClampedHessian();
  //K += _strandMesh.computeTwistingClampedHessian();
  //cout << " stretchingMine =  " << MATRIX(_strandMesh.computeStretchingClampedHessian()).format(octave) << endl;
  //cout << " bendingMine =  " << MATRIX(_strandMesh.computeBendingClampedHessian()).format(octave) << endl;
  //cout << " twistingMine =  " << MATRIX(_strandMesh.computeTwistingClampedHessian()).format(octave) << endl;
  
  VECTOR b = _dt * _dt * R;
  SPARSE_MATRIX A = _M - _dt * _dt * K;

  VECTOR rhs = _S * b;
  SPARSE_MATRIX LHS = _S * A * _S + _IminusS;

  _cgSolver.compute(LHS);
  VECTOR y = _cgSolver.solve(rhs);

  VECTOR& xDelta = y;
  _position += xDelta;
  tetStrandMesh->setDisplacement(_position.block(0, 0, _DOFsStrand, 1));
  _strandMesh.printState();
  _tetMesh.setDisplacement(_position.block(_DOFsStrand, 0, _DOFsVolume, 1));

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// indirect function so we can try out lots of different solvers
///////////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR TIMESTEPPER::solveSystem(const SPARSE_MATRIX& LHS, const VECTOR& rhs)
{
  cout << " Solving system ..." << flush;
  //VECTOR y;
  //Eigen::VectorXd y = Eigen::VectorXd::Zero(LHS.rows());
  Eigen::VectorXd y = rhs;

  bool verbose = false;
  bool usingAMG = false;
  bool usingCR = false;
  bool usingEigenCG = _disablePreconditioner;
  //bool usingEigenCG = true;

  /*
#ifdef EIGEN_USE_MKL_ALL
  {
    if (verbose)
      cout << " USING PARDISO" << endl;
    TIMER choleskyTimer("Pardiso Solve");
    Eigen::PardisoLLT<SPARSE_MATRIX> solver;  // this one is mildly faster than LDLT
    //Eigen::PardisoLDLT<SPARSE_MATRIX> solver;
    //Eigen::PardisoLU<SPARSE_MATRIX> solver;  // this one is slowest
    solver.compute(LHS);
    y = solver.solve(rhs);
    choleskyTimer.stop();
    cout << " done." << endl;
    return y;
  }

#endif
  */

  if (!_pcgEnabled)
  {
    if (verbose)
      cout << " USING CHOLESKY" << endl;
    TIMER choleskyTimer("Cholesky Solve");
#if 1
    Eigen::SimplicialLDLT<SPARSE_MATRIX> solver;
    //Eigen::SparseLU<SPARSE_MATRIX> solver;
    solver.compute(LHS);

    if (solver.info() != Eigen::Success)
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " CHOLESKY FACTORIZATION FAILED. " << endl;
      cout << " Eigenvalues: " << eigenvalues(LHS).transpose() << endl;
      exit(0);
    }
    y = solver.solve(rhs);
    if (solver.info() != Eigen::Success)
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " CHOLESKY SOLVE FAILED. " << endl;
      cout << " Eigenvalues: " << eigenvalues(LHS).transpose() << endl;
      exit(0);
    }
#else
    MATRIX A(LHS);
    y = A.colPivHouseholderQr().solve(rhs);
#endif
    choleskyTimer.stop();
    cout << " done." << endl;
    return y;
  }

  if (usingAMG)
  {
#if USING_AMGCL
    // Setup the solver:
    typedef amgcl::make_solver<
        amgcl::amg<
            amgcl::backend::eigen<double>,
            amgcl::coarsening::smoothed_aggregation,
            //amgcl::relaxation::spai0
            amgcl::relaxation::ilu0
            >,
        //amgcl::solver::bicgstab<amgcl::backend::eigen<double> >
        amgcl::solver::cg<amgcl::backend::eigen<double> >
        > Solver;

    /*
    // CG solver preconditioned with ILU0
    typedef amgcl::make_solver<
      amgcl::preconditioner::dummy<amgcl::backend::eigen<double> >,
      amgcl::solver::cg<
          amgcl::backend::eigen<double>
          >
      > Solver;
    */

    Solver::params prm;
    prm.solver.tol = 1e-6;
    prm.solver.maxiter= 1000;
    //prm.solver.verbose = true;
    prm.solver.verbose = false;

    Solver solver(LHS, prm);
    int    iters;
    double error;
    std::tie(iters, error) = solver(LHS, rhs, y);
    if (verbose)
      printf("  AMGCL iters: %3i err: %6.4e \n", iters, (float)error);
#else
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " AMG CL DISABLED. " << std::endl;
    exit(0);
#endif
  }
  else if (usingCR)
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
    cout << " Residual: " << residual.norm() << endl;
    
    //DIAGONAL diagonal(LHS);
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

// todo: solveEnergyDamped
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::solveDynamics(const bool verbose)
{
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " STRAND::TIMESTEPPER RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
    cout << " E = " << _strandMesh.E() << endl;
    cout << "Strand collisionEps = " << _strandMesh.collisionEps() << endl;
    cout << "Strand collision stiffness = " << _collisionStiffnessStrand << endl;
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
  // buildConstraintMatrix();
  buildConstraintMatrixFaster();

  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  if (tetStrandMesh == NULL)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " Wrong STRAND_MESH type passed to VOLUME_TIMESTEPPER!!!" << std::endl;
    assert(tetStrandMesh);
  }
  tetStrandMesh->computeFs();
  tetStrandMesh->computeSVDs();

  _tetMesh.setDisplacement(_position.block(_DOFsStrand, 0, _DOFsVolume, 1));
  _tetMesh.computeFs();
  _tetMesh.computeSVDs();

  //_strandMesh.setDisplacement(_position);

  // do collision detection, including spatial data structure updates 
  if (_collisionsEnabled)
    _strandMesh.computeEdgeEdgeCollisions(verbose);

  computeCollisionDetectionVolume();
  computeCollisionDetectionStrandVolume();
  // TODO: add plane constraints so that this starts firing
  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  updateConstraintTargets();
  VECTOR z =_IminusS * _constraintTargets;
  //prologueTimer.stop();

  // get the internal forces
  TIMER forceTimer("Internal forces");
  const VECTOR stretchingForce = tetStrandMesh->computeStretchingForces();
  const VECTOR bendingForce = tetStrandMesh->computeBendingForces();
  const VECTOR twistingForce = tetStrandMesh->computeTwistingForces();
  // const VECTOR stretchingForce = _strandMesh.computeStretchingForces();
  // const VECTOR bendingForce = _strandMesh.computeBendingForces();
  // const VECTOR twistingForce = _strandMesh.computeTwistingForces();
  VECTOR RVolume = _tetMesh.computeHyperelasticForces(_hyperelastic);

  if (verbose)
  {
    cout << " Bending force :   " << bendingForce.norm() << endl;
    cout << " Stretching force norm: " << stretchingForce.norm() << endl;
    // cout << " Stretching force: " << stretchingForce.block<9, 1>(0,0) << endl;
    cout << " Twisting force:   " << twistingForce.norm() << endl;
  }

  // making non-const so we can add collision forces
  VECTOR RStrand = bendingForce + stretchingForce + twistingForce;
  // VECTOR RStrand = bendingForce + twistingForce;
  VECTOR R(_DOFs);
  // cout << "R Strand norm: "<<RStrand.norm() / RStrand.rows()<<endl;
  // cout << "R Volume norm: "<<RVolume.norm() / RVolume.rows()<<endl;
  R << RStrand, RVolume;
  forceTimer.stop();

  // get the reduced stiffness matrix
  //TIMER assemblyTimer("Matrix Assembly");

  SPARSE_MATRIX KStrand;
  if (_hessianClampingEnabled)
  {
    if (verbose)
    {
      cout << " USING CLAMPED HESSIAN" << endl;
      cout << " Building material Hessian ..." << flush;
    }

    // if fast assembly is available, use it
    KStrand = tetStrandMesh->computeHyperelasticClampedHessianFast();
    if (verbose)
      cout << " done." << endl;
  }
  else
  {
    //SPARSE_MATRIX KStrand = _strandMesh.computeStretchingHessian();
    //KStrand += _strandMesh.computeBendingHessian();
    //KStrand += _strandMesh.computeTwistingHessian();
    SPARSE_MATRIX Ks = _strandMesh.computeStretchingHessian(_stretchingEnergy);
    SPARSE_MATRIX Kb = _strandMesh.computeBendingHessian();
    SPARSE_MATRIX Kt = _strandMesh.computeTwistingHessian();
    if (verbose)
    {
      cout << " Bending hessian:    " << Kb.norm() << endl;
      cout << " Stretching hessian: " << Ks.norm() << endl;
      cout << " Twisting hessian:   " << Kt.norm() << endl;
    }
    
    KStrand = Ks + Kb + Kt;
  }

  SPARSE_MATRIX KVolume = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);
  SPARSE_MATRIX K(_DOFs, _DOFs);
  concatSparseMatrix(KStrand, KVolume, K);
  cout<<"K Strand norm: "<<KStrand.norm() / KStrand.rows()/KStrand.cols()<<endl;
  cout<<"K Volume norm: "<<KVolume.norm() / KVolume.rows()/KVolume.cols()<<endl;

  //assemblyTimer.stop();

  // collision damping only appears on LHS
  const int rank = R.size();

  SPARSE_MATRIX collisionC(rank,rank);
  if (_collisionsEnabled)
  {
    // compute collision forces and stiffnesses
    // todo: strand volume case
    computeCollisionResponse(R,K,collisionC, true);
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
  
#if 1
  TIMER lhsProjectionTimer("LHS PPCG projection");
  // SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
  SPARSE_MATRIX LHS = filteredSystem();
  lhsProjectionTimer.stop();
#else
  TIMER lhsProjectionTimer("LHS PPCG projection");
  SPARSE_MATRIX AN = (_A * _N).pruned();
  SPARSE_MATRIX ANT = AN.transpose();
  
  SPARSE_MATRIX leftRight = (_N * AN).pruned();

  _A += -(AN + ANT) + leftRight + _N;   // final add takes the longest
  const SPARSE_MATRIX& LHS = _A;
  lhsProjectionTimer.stop();
#endif

  if (verbose)
  {
    cout << " b norm: " << _b.norm() << endl;
    cout << " r norm: " << c.norm() << endl;
    cout << " z norm: " << z.norm() << endl;
    cout << " targets norm: " << _constraintTargets.norm() << endl;
    cout << " rhs norm: " << rhs.norm() << endl;
    cout << " LHS norm: " << LHS.norm() << endl;
    cout << " A norm:   " << _A.norm() << endl;
    cout << " K norm:   " << K.norm() << endl;
    cout << " collisionC norm:   " << collisionC.norm() << endl;
    cout << " C norm:   " << C.norm() << endl;
    cout << " M norm:   " << _M.norm() << endl;
    cout << " IminusS norm:   " << _IminusS.norm() << endl;
  }

  // solve system using whatever solver is activated right now
  VECTOR y = solveSystem(LHS, rhs);
 
  TIMER postTimer("Epilogue");

  // aliasing _solution to \Delta x just to make clear what we're doing here
  VECTOR& xDelta = _solution;
  // TODO: add plane constraints
  xDelta = y + z;
  //xDelta = y;

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
    // buildConstraintMatrix();
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
  tetStrandMesh->setDisplacement(_position.block(0, 0, _DOFsStrand, 1));
  _tetMesh.setDisplacement(_position.block(_DOFsStrand, 0, _DOFsVolume, 1));

  // update velocity
  _velocity = invDt * (_position - _positionOld);

  // update acceleration
  _acceleration = invDt * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;
  _acceleration = _S * _acceleration;

  postTimer.stop();

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// todo: update this
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::solveNewton(const bool verbose)
{
  const REAL invDt = 1.0 / _dt;

  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " STRAND::TIMESTEPPER NEWTON SOLVE " << _currentTimestep << endl;
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

  // store the internal forces for later
  VECTOR R;
  VECTOR z;

  // do Newton-Raphson
  REAL eps = _residualTolerance;
  REAL maxR = eps * 10;
  unsigned int step = 0;

  const int minNewtonIterations = 2;

  REAL bestR = FLT_MAX;
  VECTOR bestPosition = _position;

  const REAL alpha0 = 1.0 / (_dt * _dt);
  const REAL alpha1 = 1.0 / _dt;
  const REAL alpha3 = 1.0 / _dt;
  //while (step < _maxNewtonIterations && maxR > eps)
  while ((step < _maxNewtonIterations && maxR > eps) || step < minNewtonIterations)
  {
    // build new constraints and see if we should break any
    findNewSurfaceConstraints(verbose);
    buildConstraintMatrix();

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
      _strandMesh.computeEdgeEdgeCollisions(verbose);

    // TODO: add plane constraints so that this starts firing
    // "z is a vector of the desired values for the constrained variables",
    // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
    // _constraintTargets did not project off the null directions
    updateConstraintTargets();
    z =_IminusS * _constraintTargets;

    // get the internal forces
    const VECTOR stretchingForce = tetStrandMesh->computeStretchingForces();
    const VECTOR bendingForce = tetStrandMesh->computeBendingForces();
    const VECTOR twistingForce = tetStrandMesh->computeTwistingForces();

#if 0
    cout << " stretching: " << stretchingForce.transpose() << endl;
    cout << " bending:    " << bendingForce.transpose() << endl;
    cout << " twisting:   " << twistingForce.transpose() << endl;
#endif

    // making non-const so we can add collision forces
    R = bendingForce + stretchingForce + twistingForce;
    //VECTOR R = twistingForce + stretchingForce;

    // get the reduced stiffness matrix
    TIMER assemblyTimer("Matrix Assembly");
#if 1
    SPARSE_MATRIX K = tetStrandMesh->computeHyperelasticClampedHessianFast();
    // K += _strandMesh.computeBendingClampedHessian();
    // K += _strandMesh.computeTwistingClampedHessian();
#else
    SPARSE_MATRIX K = _strandMesh.computeStretchingHessian();
    K += _strandMesh.computeBendingHessian();
    K += _strandMesh.computeTwistingHessian();
#endif
    assemblyTimer.stop();

    // collision damping only appears on LHS
    const int rank = R.size();
    SPARSE_MATRIX collisionC(rank, rank);

    // compute collision forces and stiffnesses
    computeCollisionResponse(R,K,collisionC, true);

    VECTOR acceleration = -alpha0 * (_position - _positionOld) + alpha1 * _velocityOld;
    VECTOR temp = _M * acceleration;

    // compute the RHS of the residual:
    TIMER rhsTimer("Forming Initial RHS");
    const REAL invDt = 1.0 / _dt;
    const REAL invDt2 = invDt * invDt;
    _velocity = alpha3 * (_position - _positionOld) - _velocityOld;
    _b = temp + C * _velocity + R + _externalForces;

    // collisionC does *not* appear here, since the damping only depends on
    // v^{t+1}, which only depends on \Delta x, which is the variable
    // we are solving for
    //_b = (invDt * _M - C - collisionC) * _velocity + R + _externalForces;
    rhsTimer.stop();

    // assemble system matrix A
    TIMER lhsTimer("Forming Initial LHS");
    // TODO: implement collisions
    //_A = _M * invDt2 - collisionC * invDt - K; // w/ Rayleigh
    //_A = _M * invDt2 - (C + collisionC) * invDt - K; // w/ Rayleigh
    _A = _M * alpha0 - (C + collisionC) * alpha3 - K; // w/ Rayleigh
    lhsTimer.stop();

    // in [TJM15], this is c = b - Az (page 8, top of column 2)
    TIMER projectionTimer("PPCG projection");
    // TODO: add plane constraints
    VECTOR c = _b - _A * z;

    // just for clarity, go ahead and form the RHS and LHS explicitly
    //
    // Since _S is sparse, this multipy could be accelerated significantly, 
    // but leaving it as it is for now
    VECTOR rhs = _S * c;
    SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
    projectionTimer.stop();

    maxR = rhs.squaredNorm();
    //maxR = c.squaredNorm();
    //cout << " Inf norm: " << c.lpNorm<Eigen::Infinity>() << endl;
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
    xDelta = y + z;

    // update positions
    _position += xDelta;

    // TODO: add plane constraints
    // when checking against normals, unfiltered should be negated for Newmark
    const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);

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
    tetStrandMesh->setDisplacement(_position.block(0, 0, _DOFsStrand, 1));
    cout << " Newton iteration " << step << " residual: " << maxR << endl;

    step++;
    
    if (step > minNewtonIterations && maxR < bestR)
    {
      bestR = maxR;
      bestPosition = _position;
    } 
  }

  // if this is a rotten Newton iteration, back up to a better one
  if (bestR < maxR)
  {
    _position = bestPosition;
    _strandMesh.setDisplacement(_position.block(0, 0, _DOFsStrand, 1));
    cout << " Backed out to residual " << bestR << endl;
  }

  // update velocity
  //_velocity = invDt * (_position - _positionOld);
  _velocity = alpha3 * (_position - _positionOld);

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
  TIMER functionTimer(__FUNCTION__);
  // back up current state
  _temp = _strandMesh.getDisplacement();

  // set to zero displacement
  VECTOR zero(_DOFsStrand);
  zero.setZero();
  _strandMesh.setDisplacement(zero);

  // get stiffness matrix at that state
  TET_STRAND_MESH* tetStrandMesh = dynamic_cast<TET_STRAND_MESH*>(&_strandMesh);
  tetStrandMesh->computeFs();
  tetStrandMesh->computeSVDs();
  SPARSE_MATRIX KStrand = tetStrandMesh->computeHyperelasticClampedHessianFast();

  // restore old state
  _strandMesh.setDisplacement(_temp);
  
  // build out the Rayleigh damping
  SPARSE_MATRIX CStrand = _rayleighAlphaStrand * _MStrand;
  CStrand += _rayleighBetaStrand * KStrand;
  //C += _rayleighBeta * K;

  //SPARSE_MATRIX C(_DOFs, _DOFs);
  //C.setIdentity();

  // volume part
  // back up current state

  VECTOR volumeTemp = _tetMesh.getDisplacement();

  // set to zero displacement
  VECTOR zeroVolume(_DOFsVolume);
  zeroVolume.setZero();
  _tetMesh.setDisplacement(zeroVolume);

  // get stiffness matrix at that state
  _tetMesh.computeFs();
  _tetMesh.computeSVDs();
  SPARSE_MATRIX KVolume = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);
  // restore old state
  _tetMesh.setDisplacement(volumeTemp);
  
  // build out the Rayleigh damping
  SPARSE_MATRIX CVolume = _rayleighAlphaVolume * _MVolume;
  // cout<<"CVolume size: "<<_MVolume.rows()<<" "<<_MVolume.cols() << " KVolume size: "<< KVolume.rows()<<" "<< KVolume.cols()<<endl;
  CVolume += _rayleighBetaVolume * KVolume;
  SPARSE_MATRIX C(_DOFs, _DOFs);
  concatSparseMatrix(CStrand, CVolume, C);

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

  // const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // is the vertex still inside the object? In that case, keep the constraint.
    const KINEMATIC_SHAPE* shape = constraint.shape;
    const int vertexID = constraint.vertexID;
    const VECTOR3& vertex = getMeshVertex(vertexID);
    const REAL signedDistance = shape->signedDistance(vertex);

    //bool debug = constraint.vertexID == 262;
    bool debug = false;
    if (debug)
    {
      cout << " Constraint for vertex: " << constraint.vertexID << endl;
      cout << "   Signed distance:     " << signedDistance << endl;
    }

    
    // if the distance is outside and large, move on
    if (signedDistance > 1e-6) 
    // if (signedDistance > 1e-5) 
    //if (signedDistance > 1e-4) 
    {
      constraint.isSeparating = true;
      changed = true;
      if (debug) cout << " CONSTRAINT IS OUTSIDE " << endl;
      continue;
    }
    

    // what direction is the solution pointing in?
    VECTOR3 xDirection;
    int vectorID;
    if (_strandMesh.edgeEnd())
    {
      vectorID = getPositionIndexEdgeEnd(vertexID);
    }  
    else
    {
      vectorID = getPositionIndexInterleaved(vertexID);
    }
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
    if (separationMagnitude > 1e-6)
    //if (signedDistance > 1e-5) 
    //if (signedDistance > 1e-4) 
    {
      constraint.isSeparating = true;
      changed = true;
      if (debug) cout << " SEPARATING " << endl;
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
    // if (separationMagnitude > 1e-6)   // velocity condition doesn't seem to ever get triggered
    //if (separationMagnitude > 1e-6 || velocitySeparation > 1e-6)
    //if (separationMagnitude > 1e-5 || velocitySeparation > 1e-5)
    //if (separationMagnitude > 0.0 || velocitySeparation > 0.0)
    //if (separationMagnitude > 1e-4 || velocitySeparation > 1e-4)  // jitter goes away for SNH lambda = 1000, but gets sticky
    //if (separationMagnitude > 9e-5 || velocitySeparation > 9e-5)
  //   {
  //     //cout << " Separation: " << separationMagnitude << " velocity: " << velocitySeparation << endl;
  //     constraint.isSeparating = true;
  //     changed = true;

  //     if (debug)
  //       cout << " SEPARATING" << endl;
  //   }
  }

  return changed;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// add a gravity body force to the simulation
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::addGravity(const VECTOR3& bodyForce)
{
  int totalVertices = _totalVertices;
  const vector<REAL>& oneRingVolumes = _tetMesh.restOneRingVolumes();
  if (_strandMesh.edgeEnd())
  {
    for (int x = 0; x < totalVertices; x++)
    {
      VECTOR3 scaledForce;
      if(x < _strandMesh.totalVertices())
        scaledForce = bodyForce * _strandMesh.vertexMass(x);
      else
        scaledForce = bodyForce * oneRingVolumes[x - _strandMesh.totalVertices()];
      const int index = getPositionIndexEdgeEnd(x);
      _externalForces[index]     += scaledForce[0];
      _externalForces[index + 1] += scaledForce[1];
      _externalForces[index + 2] += scaledForce[2];
    }
    return;
  }

  // const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  for (int x = 0; x < totalVertices; x++)
  {
    VECTOR3 scaledForce;
    if(x < _strandMesh.totalVertices())
      scaledForce = bodyForce * _strandMesh.vertexMass(x);
    else
      scaledForce = bodyForce * oneRingVolumes[x - _strandMesh.totalVertices()];
    const int index = getPositionIndexInterleaved(x);
    _externalForces[index]     += scaledForce[0];
    _externalForces[index + 1] += scaledForce[1];
    _externalForces[index + 2] += scaledForce[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// once attachKinematicSurfaceConstraints is called, we can also see which edges
// should be constrained
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeConstrainedEdges(const bool verbose)
{
  // hash all the vertex IDs
  map<int, bool> constrainedVertices;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++){
    if(_kinematicConstraints[x].vertexID < _strandMesh.totalVertices()){
    constrainedVertices[_kinematicConstraints[x].vertexID] = true;
    }
  }

  // go through all the edges
  _constrainedEdges.clear();
  const vector<VECTOR2I> edges = _strandMesh.edgeIndices();
  for (unsigned int x = 0; x < edges.size(); x++)
  {
    const VECTOR2I& edge = edges[x];

    // if either point is not constrained, move on
    if (constrainedVertices.find(edge[0]) == constrainedVertices.end())
      continue;
    if (constrainedVertices.find(edge[1]) == constrainedVertices.end())
      continue;

    // if both vertices are pinned, this edge should be too
    _constrainedEdges.push_back(x);
  }

  if (verbose)
    cout << " Constrained " << _constrainedEdges.size() << " edges " << endl;
}

void TIMESTEPPER::removeKinematicSurfaceConstraints(){
  _kinematicConstraints.clear();
}

void TIMESTEPPER::popBackKinematicSurfaceConstraints(){
  _kinematicConstraints.pop_back();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// constrain surface nodes inside a kinematic body to move along with that body
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::attachKinematicSurfaceConstraints(const KINEMATIC_SHAPE* shape, const bool verbose)
{
  // strand
  // get all the nodes inside the shape
  const vector<VECTOR3>& verticesStrand = _strandMesh.vertices();
  for (unsigned int x = 0; x < verticesStrand.size(); x++)
  {
    VECTOR3 v = verticesStrand[x];

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

  // computeConstrainedEdges();

  // volume
  const vector<VECTOR3>& verticesVolume = _tetMesh.vertices();
  const vector<int>& surfaceVertices = _tetMesh.surfaceVertices();
  for (unsigned int x = 0; x < surfaceVertices.size(); x++)
  {
    int whichVertex = surfaceVertices[x];
    VECTOR3 v = verticesVolume[whichVertex];

    // if it's not inside, move on
    if (!shape->inside(v)) continue;

    // if it's inside, get its local coordinates
    VECTOR3 local = shape->worldVertexToLocal(v);

    // record everything the solver will need later
    KINEMATIC_CONSTRAINT constraint;
    constraint.shape = shape;
    constraint.vertexID = whichVertex + _strandMesh.totalVertices();
    constraint.localPosition = local;

    // rememeber the constraint for later
    _kinematicConstraints.push_back(constraint);
  }

  if (verbose)
    cout << " Constrained " << _kinematicConstraints.size() << " vertices " << endl;
}

#if 1
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
    constraint.vertexID = x + _strandMesh.totalVertices();
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
  //TIMER functionTimer(__FUNCTION__);
  const vector<VECTOR3> verticesStrand = _strandMesh.vertices();
  const vector<VECTOR3> verticesVolume = _tetMesh.vertices();
  const vector<int> surfaceVerticesVolume = _tetMesh.surfaceVertices();

  if (verbose)
    cout << " Currently tracking " << _planeConstraints.size() << " plane constraints " << endl;

  // build any new constraints
  int newConstraints = 0;
  for (unsigned int y = 0; y < _collisionObjects.size(); y++)
  {
    const KINEMATIC_SHAPE* shape = _collisionObjects[y];
    // Strand
    for (unsigned int x = 0; x < verticesStrand.size(); x++)
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
      const VECTOR3& vertex = verticesStrand[vertexID];
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
    } // strand
    // volume 
    for (unsigned int x = 0; x < surfaceVerticesVolume.size(); x++)
    {
      // get the vertex
      assert(surfaceVerticesVolume[x] < int(verticesVolume.size()));
      int vertexID = surfaceVerticesVolume[x];
      int vIDGlobal = vertexID + verticesStrand.size();

      //bool debug = (vertexID == 0);
      bool debug = false;

      // if it's already in collision, skip it
      if (_inCollision[vIDGlobal]) 
      {
        if (debug) cout << " vertex is already collision, moving on" << endl;
        continue;
      }

      // see if it's inside the shape
      const VECTOR3& vertex = verticesVolume[vertexID];
      if (!shape->inside(vertex)) 
      {
        if (debug) cout << " vertex is not inside the shape, moving on" << endl;
        continue;
      }

      VECTOR3 closestPoint;
      VECTOR3 closestNormal;
      shape->getClosestPoint(vertex, closestPoint, closestNormal);
     
      // if the velocity is pulling away from the surface, don't constrain it
      VECTOR3 vertexVelocity = velocity(vIDGlobal);
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
      //if (velocitySeparation >= 0)
      //if (velocitySeparation >= -1e-9)
      if (velocitySeparation >= -FLT_EPSILON)
        continue;

      // store the constraint
      PLANE_CONSTRAINT constraint;
      constraint.shape = shape;
      constraint.vertexID = vIDGlobal;
      constraint.localClosestPoint = closestPoint;
      constraint.localNormal = closestNormal;
      constraint.isSeparating = false;
      addPlaneConstraint(constraint);

      _inCollision[vIDGlobal] = true;
      newConstraints++;
    }
  } // collision objects
  if (verbose)
    cout << " Found " << newConstraints << " new plane constraints " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the closest point positions on the surface constraints
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateSurfaceConstraints()
{
  // const vector<VECTOR3> vertices = _strandMesh.vertices();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const KINEMATIC_SHAPE& shape = *_planeConstraints[x].shape;

    // get the new vertex position
    const int vertexID = _planeConstraints[x].vertexID;
    const VECTOR3& vertex = getMeshVertex(vertexID);

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
  // if (index >= _strandMesh.totalVertices()) {
  //   VECTOR3 vertexVelocity;
  //   int idx = 3 * (index - _strandMesh.totalVertices()) + _DOFsStrand
  //   vertexVelocity[0] = _velocity[idx];
  //   vertexVelocity[1] = _velocity[idx + 1];
  //   vertexVelocity[2] = _velocity[idx + 2];
  //   return vertexVelocity;
  // }
  if (_strandMesh.edgeEnd())
  {
    assert(index >= 0);
    assert(index < _velocity.size());
    VECTOR3 vertexVelocity;
    int vectorID = getPositionIndexEdgeEnd(index);
    vertexVelocity[0] = _velocity[vectorID];
    vertexVelocity[1] = _velocity[vectorID + 1];
    vertexVelocity[2] = _velocity[vectorID + 2];
    return vertexVelocity;
  }

  // const VECTORI& globalVertexIndices = _strandMesh.globalVertexIndices();
  const int global = getPositionIndexInterleaved(index);

  VECTOR3 vertexVelocity;
  vertexVelocity[0] = _velocity[global ];
  vertexVelocity[1] = _velocity[global + 1];
  vertexVelocity[2] = _velocity[global + 2];
  return vertexVelocity;
}

#if 1
///////////////////////////////////////////////////////////////////////////////////////////////////////
// reset the Rayleigh damping constants
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::setRayelighStrand(const REAL alpha, const REAL beta)
{
  _rayleighAlphaStrand = alpha;
  _rayleighBetaStrand = beta;
}

void TIMESTEPPER::setRayelighVolume(const REAL alpha, const REAL beta)
{
  _rayleighAlphaVolume = alpha;
  _rayleighBetaVolume = beta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// do the collision detection, in anticipation of collision response
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeCollisionDetectionVolume()
{
  TIMER functionTimer(__FUNCTION__);

  // if the tet mesh has an AABB accelerator, refit it
  TET_MESH_FASTER* fast = dynamic_cast<TET_MESH_FASTER*>(&_tetMesh);
  if (fast != NULL)
    fast->refitAABB();

  // do the collision processing
  const REAL invDt = 1.0 / _dt;
  if (_vertexFaceSelfCollisionsOn)
  {
    // vertex-face collision detection
    _tetMesh.computeVertexFaceCollisions();

    // build out the vertex-face "collision tets"
    // TODO: this need to get cut down
    _tetMesh.buildVertexFaceCollisionTets(_velocity.block(_DOFsStrand, 0, _DOFsVolume, 1));
  }
  if (_edgeEdgeSelfCollisionsOn)
  {
    // edge-edge collision detection
    _tetMesh.computeEdgeEdgeCollisions();
  }
}
#endif

void TIMESTEPPER::computeCollisionDetectionStrandVolume()
{
  TIMER functionTimer(__FUNCTION__);
  // if the tet mesh has an AABB accelerator, refit it
  TET_MESH_FASTER* fast = dynamic_cast<TET_MESH_FASTER*>(&_tetMesh);
  if (fast != NULL)
    fast->refitAABB();

  // do the collision processing
  const REAL invDt = 1.0 / _dt;
  if (_vertexFaceSelfCollisionsOn)
  {
    // vertex-face collision detection  
    _tetMesh.computeVertexFaceCollisionsWithStrands(_strandMesh, _vertexFaceCollisionsStrandVolume);
    _tetMesh.buildVertexFaceCollisionTetsWithStrands(_strandMesh, _vertexFaceCollisionsStrandVolume,
                                                     _vertexFaceCollisionTetsStrandVolume,                                                   _vertexFaceCollisionAreasStrandVolume);
  }
  if (_edgeEdgeSelfCollisionsOn)
  {
    // edge-edge collision detection
    _tetMesh.computeEdgeEdgeCollisionsWithStrands(_strandMesh,
                                                  _edgeEdgeCollisionsStrandVolume, 
                                                  _edgeEdgeIntersectionsStrandVolume, 
                                                  _edgeEdgeCoordinatesStrandVolume, 
                                                  _edgeEdgeCollisionAreasStrandVolume);
  }
}

VECTOR TIMESTEPPER::computeEdgeEdgeCollisionForcesStrandVolume() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_edgeEdgeCollisionsStrandVolume.size());
  const vector<VECTOR2I>& volumeSurfaceEdges = _tetMesh.surfaceEdges();
  const vector<VECTOR2I>& strandEdges = _strandMesh.edgeIndices();
  for (unsigned int i = 0; i < _edgeEdgeCollisionsStrandVolume.size(); i++)
  {
    // Volume surface edges
    const VECTOR2I& edge0 = volumeSurfaceEdges[_edgeEdgeCollisionsStrandVolume[i].first];
    // Srtand edges
    const VECTOR2I& edge1 = strandEdges[_edgeEdgeCollisionsStrandVolume[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _tetMesh.vertices()[edge0[0]];
    vs[1] = _tetMesh.vertices()[edge0[1]];
    vs[2] = _strandMesh.vertices()[edge1[0]];
    vs[3] = _strandMesh.vertices()[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinatesStrandVolume[i].first;
    const VECTOR2& b = _edgeEdgeCoordinatesStrandVolume[i].second;

#if ADD_EDGE_EDGE_PENETRATION_BUG
    const VECTOR12 force = -_edgeEdgeCollisionAreasStrandVolume[i] * _edgeEdgeEnergyStrandVolume->gradient(vs,a,b);
#else
    const VECTOR12 force = (!_edgeEdgeIntersectionsStrandVolume[i]) ? -_edgeEdgeCollisionAreasStrandVolume[i] * _edgeEdgeEnergyStrandVolume->gradient(vs,a,b)
                                                        : -_edgeEdgeCollisionAreasStrandVolume[i] * _edgeEdgeEnergyStrandVolume->gradientNegated(vs,a,b);
#endif

    perElementForces[i] = force;
    //std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    //std::cout << " force " << i << ": " << force.norm() << " intersections: " << _edgeEdgeIntersections[i] << " area: " << _edgeEdgeCollisionAreas[i] << endl;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  VECTOR forces(_DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _edgeEdgeCollisionsStrandVolume.size(); i++)
  {
    const VECTOR2I& edge0 = volumeSurfaceEdges[_edgeEdgeCollisionsStrandVolume[i].first];
    const VECTOR2I& edge1 = strandEdges[_edgeEdgeCollisionsStrandVolume[i].second];
    const VECTOR12& edgeForce = perElementForces[i];

    vector<int> vertexIndices(4);
    // into volume
    vertexIndices[0] = edge0[0] + _strandMesh.totalVertices();
    vertexIndices[1] = edge0[1] + _strandMesh.totalVertices();
    // into strand
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int x = 0; x < 4; x++)
    {
      unsigned int index;
      if(_strandMesh.edgeEnd())
        index = getPositionIndexEdgeEnd(vertexIndices[x]);
      else
        index = getPositionIndexInterleaved(vertexIndices[x]);
      assert((int)index < _DOFs);
      forces[index]     += edgeForce[3 * x];
      forces[index + 1] += edgeForce[3 * x + 1];
      forces[index + 2] += edgeForce[3 * x + 2];
    }
  }
  
  return forces;
}

SPARSE_MATRIX TIMESTEPPER::computeEdgeEdgeCollisionClampedHessianStrandVolume() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_edgeEdgeCollisionsStrandVolume.size());
  const vector<VECTOR2I>& volumeSurfaceEdges = _tetMesh.surfaceEdges();
  const vector<VECTOR2I>& strandEdges = _strandMesh.edgeIndices();
  for (unsigned int i = 0; i < _edgeEdgeCollisionsStrandVolume.size(); i++)
  {
    // Volume surface edges
    const VECTOR2I& edge0 = volumeSurfaceEdges[_edgeEdgeCollisionsStrandVolume[i].first];
    // Srtand edges
    const VECTOR2I& edge1 = strandEdges[_edgeEdgeCollisionsStrandVolume[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _tetMesh.vertices()[edge0[0]];
    vs[1] = _tetMesh.vertices()[edge0[1]];
    vs[2] = _strandMesh.vertices()[edge1[0]];
    vs[3] = _strandMesh.vertices()[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinatesStrandVolume[i].first;
    const VECTOR2& b = _edgeEdgeCoordinatesStrandVolume[i].second;
#if ADD_EDGE_EDGE_PENETRATION_BUG
    const MATRIX12 H = -_edgeEdgeCollisionAreasStrandVolume[i] * _edgeEdgeEnergyStrandVolume->clampedHessian(vs,a,b);
#else
    const MATRIX12 H = (!_edgeEdgeIntersectionsStrandVolume[i]) ? -_edgeEdgeCollisionAreasStrandVolume[i] * _edgeEdgeEnergyStrandVolume->clampedHessian(vs,a,b)
                                                    : -_edgeEdgeCollisionAreasStrandVolume[i] * _edgeEdgeEnergyStrandVolume->clampedHessianNegated(vs,a,b);
#endif
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _edgeEdgeCollisionsStrandVolume.size(); i++)
  {
    const MATRIX12& H = perElementHessians[i];
    const VECTOR2I& edge0 = volumeSurfaceEdges[_edgeEdgeCollisionsStrandVolume[i].first];
    const VECTOR2I& edge1 = strandEdges[_edgeEdgeCollisionsStrandVolume[i].second];

    vector<int> vertexIndices(4);
    // into volume
    vertexIndices[0] = edge0[0] + _strandMesh.totalVertices();
    vertexIndices[1] = edge0[1] + _strandMesh.totalVertices();
    // into strand
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int y = 0; y < 4; y++)
    {
      int yVertex = vertexIndices[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = vertexIndices[x];
        unsigned int indexX, indexY;
        if(_strandMesh.edgeEnd()){
          indexX = getPositionIndexEdgeEnd(xVertex);
          indexY = getPositionIndexEdgeEnd(yVertex);
        }
        else{
          indexX = getPositionIndexInterleaved(xVertex);
          indexY = getPositionIndexInterleaved(yVertex);
        }
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(indexX + a, indexY + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }
 
  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

VECTOR TIMESTEPPER::computeVertexFaceCollisionForcesStrandVolume() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_vertexFaceCollisionTetsStrandVolume.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTetsStrandVolume.size(); i++)
  {
    vector<VECTOR3> vs(4);
    // on strand
    vs[0] = _strandMesh.vertices()[_vertexFaceCollisionTetsStrandVolume[i][0]];
    for (unsigned int j = 1; j < 4; j++)
      vs[j] = _tetMesh.vertices()[_vertexFaceCollisionTetsStrandVolume[i][j]];
    const VECTOR12 force = -_vertexFaceCollisionAreasStrandVolume[i] * _vertexFaceEnergyStrandVolume->gradient(vs);
    perElementForces[i] = force;

#if ENABLE_DEBUG_TRAPS
    if (force.hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " NaN in collision tet: " << i << endl;
      for (int j = 0; j < 4; j++)
        cout << " v" << j << ": " << vs[j].transpose() << endl;
      cout << " gradient: " << endl << _vertexFaceEnergyStrandVolume->gradient(vs) << endl;
    }
#endif
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  VECTOR forces(_DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _vertexFaceCollisionTetsStrandVolume.size(); i++)
  {
    VECTOR4I tet = _vertexFaceCollisionTetsStrandVolume[i];
    for (int i = 1; i < 4; i++)
    {
      tet[i] += _strandMesh.totalVertices();
    }
    const VECTOR12& tetForce = perElementForces[i];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index;
      if(_strandMesh.edgeEnd())
        index = getPositionIndexEdgeEnd(tet[x]);
      else
        index = getPositionIndexInterleaved(tet[x]);
      // unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

SPARSE_MATRIX TIMESTEPPER::computeVertexFaceCollisionClampedHessianStrandVolume() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_vertexFaceCollisionTetsStrandVolume.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTetsStrandVolume.size(); i++)
  {
    vector<VECTOR3> vs(4);
    // on strand
    vs[0] = _strandMesh.vertices()[_vertexFaceCollisionTetsStrandVolume[i][0]];
    for (unsigned int j = 1; j < 4; j++)
      vs[j] = _tetMesh.vertices()[_vertexFaceCollisionTetsStrandVolume[i][j]];
    const MATRIX12 H = -_vertexFaceCollisionAreasStrandVolume[i] * _vertexFaceEnergyStrandVolume->clampedHessian(vs);
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _vertexFaceCollisionTetsStrandVolume.size(); i++)
  {
    VECTOR4I tet = _vertexFaceCollisionTetsStrandVolume[i];
    for (int i = 1; i < 4; i++)
    {
      tet[i] += _strandMesh.totalVertices();
    }
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        unsigned int indexX, indexY;
        if(_strandMesh.edgeEnd()){
          indexX = getPositionIndexEdgeEnd(xVertex);
          indexY = getPositionIndexEdgeEnd(yVertex);
        }
        else{
          indexX = getPositionIndexInterleaved(xVertex);
          indexY = getPositionIndexInterleaved(yVertex);
        }
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(indexX + a, indexY + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
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
  VECTOR collisionForcesStrand(_DOFsStrand);
  SPARSE_MATRIX collisionKStrand(_DOFsStrand, _DOFsStrand);
  SPARSE_MATRIX collisionCStrand(_DOFsStrand, _DOFsStrand);
  VECTOR collisionForcesVolume(_DOFsVolume);
  SPARSE_MATRIX collisionKVolume(_DOFsVolume, _DOFsVolume);
  SPARSE_MATRIX collisionCVolume(_DOFsVolume, _DOFsVolume);
  VECTOR collisionForcesStrandVolume(_DOFs);
  SPARSE_MATRIX collisionKStrandVolume(_DOFs, _DOFs);
  SPARSE_MATRIX collisionCStrandVolume(_DOFs, _DOFs);

  collisionForces.setZero();
  collisionK.setZero();
  collisionC.setZero();
  collisionForcesStrand.setZero();
  collisionKStrand.setZero();
  collisionCStrand.setZero();
  collisionForcesVolume.setZero();
  collisionKVolume.setZero();
  collisionCVolume.setZero();
  collisionForcesStrandVolume.setZero();
  collisionKStrandVolume.setZero();
  collisionCStrandVolume.setZero();

  const REAL dampingBetaStrand = _collisionDampingBetaStrand;
  _strandMesh.setCollisionStiffness(_collisionStiffnessStrand);
  const REAL dampingBetaVolume = _collisionDampingBetaVolume;
  _tetMesh.setCollisionStiffness(_collisionStiffnessVolume);
  const REAL dampingBetaStrandVolume = _collisionDampingBetaStrandVolume;
  _vertexFaceEnergyStrandVolume -> mu() = _collisionStiffnessStrandVolume;
  _edgeEdgeEnergyStrandVolume -> mu() = _collisionStiffnessStrandVolume;

  //Build the collision penalties
  vector<REAL> springArgsStrand;
  springArgsStrand.push_back(_collisionStiffnessStrand);
  springArgsStrand.push_back(_collisionEps);
  ENERGY_1D* normalSpringStrand = new ENERGY_1D(springArgsStrand);
  SIGNED_LEN_PLANES* lengthFuncStrand = new SIGNED_LEN_PLANES();
  C_PLANES* cFuncStrand = new C_PLANES();

  SIGNED_LEN_PLANES* lengthFuncEEStrand = new SIGNED_LEN_PLANES();
  ENERGY_1D* normalSpringEEStrand = new ENERGY_1D(springArgsStrand);
  C_PLANES_EE* cFuncEEStrand = new C_PLANES_EE();

  ENERGY_12D* vfGeneralStrand = new ENERGY_12D(normalSpringStrand, cFuncStrand, lengthFuncStrand);
  ENERGY_12D* eeGeneralStrand = new ENERGY_12D(normalSpringEEStrand, cFuncEEStrand, lengthFuncEEStrand);

  vector<REAL> springArgsVolume;
  springArgsVolume.push_back(_collisionStiffnessVolume);
  springArgsVolume.push_back(_collisionEps);
  ENERGY_1D* normalSpringVolume = new ENERGY_1D(springArgsVolume);
  SIGNED_LEN_PLANES* lengthFuncVolume = new SIGNED_LEN_PLANES();
  C_PLANES* cFuncVolume = new C_PLANES();

  SIGNED_LEN_PLANES* lengthFuncEEVolume = new SIGNED_LEN_PLANES();
  ENERGY_1D* normalSpringEEVolume = new ENERGY_1D(springArgsVolume);
  C_PLANES_EE* cFuncEEVolume = new C_PLANES_EE();

  ENERGY_12D* vfGeneralVolume = new ENERGY_12D(normalSpringVolume, cFuncVolume, lengthFuncVolume);
  ENERGY_12D* eeGeneralVolume = new ENERGY_12D(normalSpringEEVolume, cFuncEEVolume, lengthFuncEEVolume);


  // edge-edge case strand
  VECTOR forcesEEStrand;
  SPARSE_MATRIX hessianEEStrand;
  
  forcesEEStrand = _strandMesh.computeEdgeEdgeCollisionForces(*eeGeneralStrand);
  hessianEEStrand = _strandMesh.computeEdgeEdgeCollisionClampedHessian(*eeGeneralStrand);
  
  collisionForcesStrand += forcesEEStrand;
  collisionKStrand += hessianEEStrand;

  if (dampingBetaStrand > 0.0)
    collisionCStrand += dampingBetaStrand * hessianEEStrand;

  // vertex-face volume
  VECTOR forcesVFVolume;
  SPARSE_MATRIX hessianVFVolume;
  if (_vertexFaceSelfCollisionsOn)
  {
    // get vertex-face collision forces and gradient
    forcesVFVolume = _tetMesh.computeVertexFaceCollisionForces(*vfGeneralStrand);
    hessianVFVolume = _tetMesh.computeVertexFaceCollisionClampedHessian(*vfGeneralStrand);
    
    collisionForcesVolume += forcesVFVolume;
    collisionKVolume += hessianVFVolume;
    collisionCVolume += dampingBetaVolume * hessianVFVolume;
  }

  // edge-edge volume
  VECTOR forcesEEVolume;
  SPARSE_MATRIX hessianEEVolume;
  if (_edgeEdgeSelfCollisionsOn)
  {
    forcesEEVolume = _tetMesh.computeEdgeEdgeCollisionForces(*eeGeneralStrand);
    hessianEEVolume = _tetMesh.computeEdgeEdgeCollisionClampedHessian(*eeGeneralStrand);
    
    collisionForcesVolume += forcesEEVolume;
    collisionKVolume += hessianEEVolume;
    collisionCVolume += dampingBetaVolume * hessianEEVolume;
  }

  // vertex-face strand volume
  VECTOR forcesVFStrandVolume;
  SPARSE_MATRIX hessianVFStrandVolume;
  if (_vertexFaceSelfCollisionsOn)
  {
    // get vertex-face collision forces and gradient
    forcesVFStrandVolume = computeVertexFaceCollisionForcesStrandVolume();
    hessianVFStrandVolume = computeVertexFaceCollisionClampedHessianStrandVolume();
    
    collisionForcesStrandVolume += forcesVFStrandVolume;
    collisionKStrandVolume += hessianVFStrandVolume;
    collisionCStrandVolume += dampingBetaStrandVolume * hessianVFStrandVolume;
  }

  // edge-edge strand volume
  VECTOR forcesEEStrandVolume;
  SPARSE_MATRIX hessianEEStrandVolume;
  if (_edgeEdgeSelfCollisionsOn)
  {
    forcesEEStrandVolume = computeEdgeEdgeCollisionForcesStrandVolume();
    hessianEEStrandVolume = computeEdgeEdgeCollisionClampedHessianStrandVolume();
    
    collisionForcesStrandVolume += forcesEEStrandVolume;
    collisionKStrandVolume += hessianEEStrandVolume;
    collisionCStrandVolume += dampingBetaStrandVolume * hessianEEStrandVolume;
  }

  collisionForces << collisionForcesStrand, collisionForcesVolume;
  cout << "Strand Self collision forces: " << collisionForcesStrand.norm() << endl;
  collisionForces += collisionForcesStrandVolume;
  cout << "Added strand volume collision forces: " << collisionForces.norm() << endl;
  concatSparseMatrix(collisionKStrand, collisionKVolume, collisionK);
  collisionK += collisionKStrandVolume;
  concatSparseMatrix(collisionCStrand, collisionCVolume, collisionC);
  collisionC += collisionCStrandVolume;

#if VERY_VERBOSE
  if (verbose)
  {
    cout << " collision forces: " << collisionForces.norm() << endl;
    cout << " collision K: " << collisionK.norm() << endl;
    cout << " collision C: " << collisionC.norm() << endl;
  }
#endif

  R += collisionForces;
  K += collisionK;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// apply a rigid rotation as a warm start to the solve
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::applyRigidRotation(const MATRIX3& R, const VECTOR3& translation)
{
  // see which vertices are already constrained, because they'll
  // automatically be rotated
  vector<VECTOR3>& vertices = _strandMesh.vertices();
  /*
  vector<bool> isConstrained(vertices.size());
  for (unsigned int x = 0; x < vertices.size(); x++)
    isConstrained[x] = false;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];
    const int vertexID = constraint.vertexID;
    isConstrained[vertexID] = true;
  }
  */

  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // why does uncommenting this break the warm start?
    //if (isConstrained[x]) continue;

    vertices[x] = R * vertices[x] + translation;
  }

  // update integrator quantities
  _position.block(0, 0, _DOFsStrand, 1) = _strandMesh.getDisplacement();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// see if a specific vertex is kinematically constrained
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::isKinematicallyConstrained(const int vertexID)
{
  assert(vertexID >= 0);
  assert(vertexID < _totalVertices);

  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];
    if (vertexID == constraint.vertexID)
      return true;
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildPlaneConstraintsOnly()
{
  TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsP;

  vector<bool> diagonalSeen(_DOFs);
  for (int x = 0; x < _DOFs; x++)
    diagonalSeen[x] = false;

  int totalVertices = _totalVertices;
  vector<bool> isKinematic(totalVertices);
  for (int x = 0; x < totalVertices; x++)
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
SPARSE_MATRIX TIMESTEPPER::buildPlaneConstraintsNoIdentity()
{
  // timing is negligible
  //TIMER functionTimer(__FUNCTION__);
  
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> tripletsP;

  int totalVertices = _totalVertices;
  vector<bool> isKinematic(totalVertices);
  for (int x = 0; x < totalVertices; x++)
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
// apply PPCG filters to system matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& TIMESTEPPER::filteredSystem()
{
  TIMER functionTimer("PPCG: filteredSystem");
#if 1

  // apply plane contraints
  _P = buildPlaneConstraintsOnly();
  _PnoI = buildPlaneConstraintsNoIdentity();

  // see which kinematic entries are zero
  int totalVertices = _totalVertices;
  vector<bool> isKinematic(3 * totalVertices);
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
  vector<bool> isPlanar(3 * totalVertices);
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
} // STRAND_VOLUME
