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
#ifndef TET_WISP_MESH_H
#define TET_WISP_MESH_H

#include "TET_STRAND_MESH.h"

namespace HOBAK {

using namespace std;

class TET_WISP_MESH : public TET_STRAND_MESH
{
public:
  TET_WISP_MESH() = default;

  // accepts a vector of individual strands
  TET_WISP_MESH(const vector<VECTOR3>& restVertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& baseRadius,
              const REAL& tipRadius);

  virtual ~TET_WISP_MESH();
  const VECTOR& vertexRadii() const { return _vertexRadii; };
  const VECTOR& strandLengths() const { return _strandLengths; };

protected:
  void initializeWisps();
  void computeWispMasses();
  void computeStrandIDs();

  const REAL computeTetVolume(const int index) const;

  VECTOR _vertexRadii;   // _totalVertices
  VECTOR _edgeRadii;     // _totalEdges
  VECTOR _strandLengths; // _totalStrands

  VECTORI _edgeStrandID; // _totalEdges

  REAL _baseRadius;
  REAL _tipRadius;
};

}

#endif
