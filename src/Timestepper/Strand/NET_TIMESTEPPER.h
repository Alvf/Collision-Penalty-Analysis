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
#ifndef STRAND_VOLUME_TIMESTEPPER_H
#define STRAND_VOLUME_TIMESTEPPER_H

#include "VOLUME_TIMESTEPPER.h"
#include "Geometry/TET_STRAND_MESH.h"

#include <iostream>

namespace HOBAK {
namespace STRAND {

////////////////////////////////////////////////////////////////////////////////////////////////////
// This implements Baraff-Witkin-style constraints from
//  "Large Steps in Cloth Simulation", SIGGRAPH 1998
// by building the system described in 
//  "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class NET_TIMESTEPPER : public VOLUME_TIMESTEPPER
{
public:
  NET_TIMESTEPPER(STRAND_NET_MESH& strandMesh, STRAND::STRETCHING& stretching);
  virtual ~NET_TIMESTEPPER();

  // take a timestep
  virtual bool solveDynamics(const bool verbose) override;
  virtual bool solveDynamicsWithRotation(const bool verbose, const MATRIX3& rotation, const VECTOR3& translation);

  void stiffenKinematicStrands();

protected:
  virtual void buildConstraintMatrix() override;
  virtual void buildConstraintMatrixFaster() override;
  
  // build the damping matrix based on the rest pose stiffness
  virtual SPARSE_MATRIX buildRayleighDampingMatrix() override;
  
  // indirect function so we can try out lots of different solvers
  virtual VECTOR solveSystem(const SPARSE_MATRIX& LHS, const VECTOR& rhs);

  // apply PPCG filters to system matrix
  SPARSE_MATRIX& filteredSystem();

  // build plane-only constraint matrix
  SPARSE_MATRIX buildPlaneConstraintsOnly(); 
  SPARSE_MATRIX buildPlaneConstraintsNoIdentity();

  // preallocated sparse matrices for fast constraint projection 
  SPARSE_MATRIX _BA, _C, _P, _PnoI, _PAP;

  // preallocated sparse matrices for missing term in fast constraint projection 
  SPARSE_MATRIX _missing, _blockDiag;
};

} // HOBAK
} // NET_TIMESTEPPER

#endif
