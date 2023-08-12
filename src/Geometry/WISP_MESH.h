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
#ifndef WISP_MESH_H
#define WISP_MESH_H

#include "STRAND_MESH.h"

namespace HOBAK {

using namespace std;

class WISP_MESH : public STRAND_MESH
{
public:
  WISP_MESH() = default;

  // accepts a vector of individual strands
  WISP_MESH(const vector<VECTOR3>& restVertices,
            const vector<vector<int> >& strandIndices,
            const REAL& E,
            const REAL& G,
            const REAL& density,
            const REAL& baseRadius,
            const REAL& tipRadius);

  virtual ~WISP_MESH();
  const VECTOR& wispRadii() const { return _wispRadii; };
  const VECTOR& strandLengths() const { return _strandLengths; };

  // compute stretching quantities
  virtual REAL computeStretchingEnergy(const STRAND::STRETCHING& hyperelastic) const override;
  virtual VECTOR computeStretchingForces(const STRAND::STRETCHING& hyperelastic) const override;
  virtual SPARSE_MATRIX computeStretchingHessian(const STRAND::STRETCHING& hyperelastic) const override;
  virtual SPARSE_MATRIX computeStretchingClampedHessian(const STRAND::STRETCHING& hyperelastic) const override;

protected:
  // generic initialization across multiple constructors
  void initializeWisps();

  // compute variable-density wisp masses
  void computeWispMasses();

  VECTOR _wispRadii; // _totalVertices

  VECTOR _strandLengths; // _totalStrands

  VECTOR _kStretch; // _totalEdges

  // different radii along the strand
  REAL _baseRadius;
  REAL _tipRadius;
};

}

#endif
