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
#ifndef TRIANGLE_MESH_FASTER_H
#define TRIANGLE_MESH_FASTER_H

#include "TRIANGLE_MESH.h"
#include "AABB_TREE.h"

namespace HOBAK {

using namespace std;

class TRIANGLE_MESH_FASTER: public TRIANGLE_MESH
{
public:
  TRIANGLE_MESH_FASTER() = default;
  TRIANGLE_MESH_FASTER(const vector<VECTOR3>& restVertices, 
                const vector<VECTOR3I>& triangles);

  // find all the vertex-face collision pairs, using the InFaceRegion test
  virtual void computeVertexFaceCollisions() override;
  
  // find all the edge-edge collision pairs
  virtual void computeEdgeEdgeCollisions() override;

protected:
  // mapping from edge index pairs to _surfaceEdges
  map<pair<int, int>, int> _edgeHash;

  // collision detection acceleration structure for triangles
  AABB_TREE _aabbTreeTriangles;
  
  // collision detection acceleration structure for edges
  AABB_TREE _aabbTreeEdges;
};

}

#endif
