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
#ifndef STRAND_MESH_FASTER_H
#define STRAND_MESH_FASTER_H

#include "STRAND_MESH.h"
#include "Geometry/AABB_TREE.h"

namespace HOBAK {

using namespace std;

class STRAND_MESH_FASTER : public STRAND_MESH
{
public:
  // accepts a vector of individual strands
  STRAND_MESH_FASTER(const vector<VECTOR3>& restVertices,
                  const vector<vector<int> >& strandIndices,
                  const REAL& E,        // Young's modulus
                  const REAL& nu,       // Poissons' ratio
                  const REAL& density,
                  const REAL& radiusA,
                  const REAL& radiusB);

  virtual ~STRAND_MESH_FASTER();

  // collision detection
  virtual void computeEdgeEdgeCollisions(const REAL& collisionEps, const bool verbose) override;

  // twisting energy
  virtual SPARSE_MATRIX computeTanTwistingHessian() const override;
  virtual SPARSE_MATRIX computeTanTwistingClampedHessian() const override;

  // compute stretching quantities
  virtual SPARSE_MATRIX computeStretchingHessian(const STRAND::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingClampedHessian(const STRAND::STRETCHING& hyperelastic) const;

  // fast computation of all the clamped Hessians for the elasticity energies
  SPARSE_MATRIX computeClampedElasticityHessian(const STRAND::STRETCHING& stretching);

  // compute collision quantities
  virtual SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessian(const ENERGY_12D& edgeEdgeEnergy) const override;

protected:
  // bending energy
  virtual SPARSE_MATRIX computeTanBendingHessian();
  virtual SPARSE_MATRIX computeTanBendingClampedHessian();

  // bake out the sparsity pattern
  void computeMatrixSparsity();
  void computeMatrixSparsityEdgeEnd();
  void computeMatrixSparsityInterleaved();

  // find the compressed index mapping
  void computeCompressedIndices();
  void computeCompressedIndicesEdgeEnd();
  void computeCompressedIndicesInterleaved();

  // shared matrix construction between clamped and unclamped
  //SPARSE_MATRIX buildPerBendMatrix(const vector<MATRIX11>& perBendHessians) const;

  // for sparse matrix entry (x,y), find the compressed index
  map<pair<int,int>, int> _compressedIndex;
  
  // global sparsity pattern
  SPARSE_MATRIX _sparseK;

  // for each entry in the global stiffness matrix, the
  // tet indices to gather entries from
  vector<vector<VECTOR3I> > _hessianBendGathers;
  vector<vector<VECTOR3I> > _hessianEdgeGathers;

  // cache the per-bend hessians
  vector<MATRIX11> _perBendHessians;
  
  // cache the per-bend hessians
  vector<MATRIX6> _perEdgeHessians;

  AABB_TREE* _collisionTree;

};

}

#endif
