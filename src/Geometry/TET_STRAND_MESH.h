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
#ifndef TET_STRAND_MESH_H
#define TET_STRAND_MESH_H

#include "STRAND_MESH.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Hyperelastic/Strand/COMPOSITE.h"
#include "Geometry/AABB_TREE.h"
#include <array>

namespace HOBAK {

using namespace std;

class TET_STRAND_MESH : public STRAND_MESH
{
public:
  TET_STRAND_MESH() = default;

  // accepts a vector of individual strands
  TET_STRAND_MESH(const vector<VECTOR3>& restVertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  // accepts a vector of individual strands, and a modified initial pose
  TET_STRAND_MESH(const vector<VECTOR3>& restVertices,
              const vector<VECTOR3>& vertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  vector<VECTOR4I>& tets()                  { return _tets; };
  const vector<VECTOR4I>& tets() const      { return _tets; };
  vector<VECTOR3I>& tetEdges()              { return _tetEdges; };
  const vector<VECTOR3I>& tetEdges() const  { return _tetEdges; };
  const int totalTets() const               { return _totalTets; };
  virtual const int DOFs() const override   { return _vertices.size() * 3; };
  const vector<bool>& isTetInverted() const { return _isTetInverted; };
  const vector<MATRIX3>& volumeFs() const { return _volumeFs; };
  vector<STRAND::COMPOSITE>& materials()    { return _materials; };
  const vector<STRAND::COMPOSITE>& materials() const { return _materials; };
  const vector<VECTOR12>& perTetForces() const { return _perTetForces; };
  bool& strandSelfCollisionDisabled() { return _strandSelfCollisionDisabled; };

  // scale strand material parameters
  void scaleTwistingMu(const REAL& scale);
  void scaleBendingMu(const REAL& scale);
  void scaleStretchingMu(const REAL& scale);

  virtual ~TET_STRAND_MESH();

  void computeFs();
  void computeSVDs();

  VECTOR computeHyperelasticForces();
  SPARSE_MATRIX computeHyperelasticClampedHessian();
  SPARSE_MATRIX computeHyperelasticClampedHessianFast();

  VECTOR computeStretchingForces();
  VECTOR computeBendingForces();
  VECTOR computeTwistingForces();
  
  REAL computeStretchingEnergy();
  REAL computeBendingEnergy();
  REAL computeTwistingEnergy();

  virtual void setDisplacement(const VECTOR& delta) override;
  virtual const VECTOR getDisplacement() const override;

  // collision detection
  virtual void computeEdgeEdgeCollisions(const REAL& collisionEps, const bool verbose) override;

protected:
  void initializeTets();
  void initializeFastHessian();
  void computeTetIndices();
  void computeDmInvs();
  void computePFpxs();
  void computeTwistPFpxs();

  // compute the volume of a tet
  static REAL computeTetVolume(const vector<VECTOR3>& tetVertices);

  // compute volumes for tets -- works for rest and deformed, just
  // pass it _restVertices or _vertices
  VECTOR computeTetVolumes(const vector<VECTOR3>& vertices);

  void recomputeTwistNormals();
  
  // find the compressed index mapping
  void computeCompressedIndices();

  // compute the vertex and edge indices inside the global vector
  void rebuildGlobalIndices();

  // overload the collision matrix
  virtual SPARSE_MATRIX buildEdgeEdgeMatrix(const vector<MATRIX12>& perEdgeHessians) const override;

  // find which strand each edge belongs to
  void buildEdgeStrandIndices();

  // tet-specific quantities
  vector<VECTOR4I> _tets;
  vector<VECTOR3I> _tetEdges; // which edges are contained by this tet?
  vector<VECTOR2I> _tetBends; // which bends are contained by this tet?
  VECTOR  _restTetVolumes;
  int _totalTets;
  
  // support for computing volume deformation gradient F
  vector<MATRIX3> _volumeDmInvs;
  vector<MATRIX3> _Us;
  vector<VECTOR3> _Sigmas;
  vector<MATRIX3> _Vs;

  // support for computing edge deformation gradients
  VECTOR _edgeDmInvs;
  
  // deformation gradients
  vector<MATRIX3> _volumeFs;
  vector<MATRIX3> _twistFs; //Fs for twist energy
  vector<pair<REAL,REAL>> _abWeights; //alpha/beta for twist energy
  vector<VECTOR3> _edgeFs;
  vector<MATRIX3x2> _bendingEs;

  // change-of-basis to go from deformation gradient (F) to positions (x)
  vector<MATRIX9x12> _volumePFPxs;
  vector<MATRIX9x12> _twistPFPxs;

  vector<STRAND::COMPOSITE> _materials;

  // is the vertex at the beginning or end of a strand?
  vector<bool> _isStrandEnd;   // _totalVertices
  vector<bool> _isStrandBegin; // _totalVertices

  // which tets are inverted?
  vector<bool> _isTetInverted; // _totalTets

  // fast Hessian construction
  mutable SPARSE_MATRIX _sparseA;

  // for sparse matrix entry (x,y), find the compressed index
  map<pair<int,int>, int> _compressedIndex;
  
  // for each entry in the global stiffness matrix, the
  // tet indices to gather entries from
  vector<vector<VECTOR3I> > _hessianGathers;

  AABB_TREE* _collisionTree;

  // DEBUG: peek at forces
  vector<VECTOR12> _perTetForces;

  // are intra-strand collisions allowed?
  bool _strandSelfCollisionDisabled;

  // strand index of each edge
  VECTOR _perEdgeStrandIndex;
};

}

#endif
