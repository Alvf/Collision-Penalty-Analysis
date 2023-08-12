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
#ifndef STRAND_NET_MESH_H
#define STRAND_NET_MESH_H

#include "STRAND_MESH.h"

namespace HOBAK {

using namespace std;

class STRAND_NET_MESH : public STRAND_MESH
{
public:
  STRAND_NET_MESH() = default;

  // accepts a vector of individual strands
  STRAND_NET_MESH(const vector<VECTOR3>& restVertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);
  STRAND_NET_MESH(const vector<VECTOR3>& restVertices,
              const vector<VECTOR3>& vertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  virtual ~STRAND_NET_MESH();

  // stretching energy
  virtual VECTOR computeStretchingForces(const STRAND::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingHessian(const STRAND::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingClampedHessian(const STRAND::STRETCHING& hyperelastic) const;

  // bending energy
  VECTOR computeBendingForces();
  SPARSE_MATRIX computeBendingHessian();
  SPARSE_MATRIX computeBendingClampedHessian();

  // twisting energy
  VECTOR computeTwistingForces() const;
  SPARSE_MATRIX computeTwistingHessian() const;
  SPARSE_MATRIX computeTwistingClampedHessian() const;

  virtual void setDisplacement(const VECTOR& delta);
  virtual const VECTOR getDisplacement() const;
  void updateProperties();
  virtual const int DOFs() const    { return _vertices.size() * 3;};

// too jittery for now
#if 0
  // collisions that take into account wisp radii
  virtual VECTOR computeEdgeEdgeCollisionForces() const override;
  virtual SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessian() const override;
  virtual void computeEdgeEdgeCollisions(const bool verbose = false) override;
#endif

protected:
  virtual void initialize() override;
  REAL clamp(const REAL& input, const REAL& bottom, const REAL& top);
  void buildGlobalIndices();
  VECTOR buildPerBendVector(const vector<VECTOR9>& perBendForces) const;
  SPARSE_MATRIX buildPerBendMatrix(const vector<MATRIX9>& perBendHessians) const;
  VECTOR buildPerEdgeVector(const vector<VECTOR6>& perEdgeForces) const;
  SPARSE_MATRIX buildPerEdgeMatrix(const vector<MATRIX6>& perEdgeHessians) const;
  VECTOR3 findOrthogonal(const VECTOR3& u);
  // vertex to (the other vertex on the edge, edgeIdx)
  vector<vector<VECTOR2I>> _vertexToEdge;
  VECTOR _restThetas; 
  MATRIX6x9 _pFpx;
  // MATRIX6x9 _edgepFpx;

  vector<bool> _isEnd;

};

}

#endif
