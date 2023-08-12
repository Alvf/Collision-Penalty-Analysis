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
#ifndef TET_MESH_FASTER_H
#define TET_MESH_FASTER_H

#include "TET_MESH.h"
#include "AABB_TREE.h"

namespace HOBAK {

class TET_MESH_FASTER : public TET_MESH
{
public:
  TET_MESH_FASTER(const vector<VECTOR3>& restVertices, 
                  const vector<VECTOR4I>& tets);

  // do something so that this doesn't run so slow
  virtual SPARSE_MATRIX computeHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const override;
  virtual SPARSE_MATRIX computeDampingHessian(const VOLUME::DAMPING& damping) const override;
  virtual BLOCK_SPARSE_MATRIX3 computeBlockHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const;

  // find all the vertex-face collision pairs, using the InFaceRegion test
  virtual void computeVertexFaceCollisions(const REAL& collisionEps) override;

  // find all the edge-edge collision pairs
  virtual void computeEdgeEdgeCollisions(const REAL& collisionEps) override;

  const AABB_TREE& aabbTreeTriangles() const { return _aabbTreeTriangles; };
  void refitAABB() { _aabbTreeTriangles.refit(); _aabbTreeEdges.refit(); };

private:
  // find the compressed index mapping
  void computeCompressedIndices();

  mutable bool _sparsityCached;
  mutable SPARSE_MATRIX _sparseA;

  mutable bool _blockSparsityCached;
  //mutable BLOCK_SPARSE_MATRIX3 _blockSparseA;

  // for sparse matrix entry (x,y), find the compressed index
  map<pair<int,int>, int> _compressedIndex;

  // cache the hessian for each tet
  mutable vector<MATRIX12> _perElementHessians;

  // mapping from edge index pairs to _surfaceEdges
  map<pair<int, int>, int> _edgeHash;

  // for each entry in the global stiffness matrix, the
  // tet indices to gather entries from
  vector<vector<VECTOR3I> > _hessianGathers;

  // collision detection acceleration structure for triangles
  AABB_TREE _aabbTreeTriangles;
  
  // collision detection acceleration structure for edges
  AABB_TREE _aabbTreeEdges;
};

}

#endif
