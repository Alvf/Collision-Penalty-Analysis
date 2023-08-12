#ifndef FORK_CCD_TEST_H
#define FORK_CCD_TEST_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

class FORK_CCD_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " A multi-pronged mesh is dropped on the ground, and the collisions " << endl;
    cout << " are sufficiently severe that early on there is flickering because " << endl;
    cout << " the penetration goes past the collision epsilon in a single timetep" << endl;
    cout << endl;
    cout << " CCD is (presumably) needed for this scenario, but NOT IMPLEMENTED." << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "fork_ccd";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/fork/fork_25.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);
    vertices = TET_MESH::normalizeVertices(vertices);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }
    _gravity = VECTOR3(0, -1.0, 0);

    REAL E = 10.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::ARAP(mu, 1.0);

    const vector<REAL>& areas = _tetMesh->surfaceTriangleAreas();
    REAL smallest = areas[0];
    REAL largest = areas[0];
    for (unsigned int x = 1; x < areas.size(); x++)
    {
      if (areas[x] > largest) largest  = areas[x];
      if (areas[x] < largest) smallest = areas[x];
    }
    cout << "Largest triangle area:  "  << largest << endl;
    cout << "Smallest triangle area: "  << smallest << endl;

    setCollisions();

    // build the time integrator
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);

    VECTOR3 center0(0.0,-4.77, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    // collision constants
    const REAL collisionMu = 1000.0;
    _volumeSolver->collisionStiffness() = collisionMu;
    _volumeSolver->collisionDampingBeta() = 0.01;

    _volumeSolver->vertexFaceSelfCollisionsOn() = true;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = true;

    // for the crumpling E
    _eye    = VECTOR3(0.929159, 0.720205, 3.59228);
    _lookAt = VECTOR3(0.870294, 0.670356, 2.59526);
    _up     = VECTOR3(0.00648923, 0.99871, -0.0503162);
    return true;
  }
  
};

}

#endif
