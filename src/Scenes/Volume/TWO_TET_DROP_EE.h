#ifndef TWO_TET_DROP_EE_H
#define TWO_TET_DROP_EE_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

class TWO_TET_DROP_EE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " One tet is dropped on another in an edge-edge collision scenario." << endl;
    cout << " If the collision processing is working correctly, the top tet will " << endl;
    cout << " just balance on top of the bottom tet." << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "two_tet_drop_ee";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/two_tets_edge_edge.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }
    _gravity = VECTOR3(0, -0.3, 0);

    REAL E = 1.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(mu, lambda);

    // setup collision objects
    setCollisions(10.0, 0.01);

    // build the time integrator
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);
    _volumeSolver->setDt(1.0 / 30.0);

    VECTOR3 center0(0.0,-4.90, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    // collision constants
    //const REAL collisionMu = 10.0;
    //_volumeSolver->collisionStiffness() = collisionMu;
    _volumeSolver->collisionDampingBeta() = 0.01;

    _volumeSolver->vertexFaceSelfCollisionsOn() = false;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = true;

    _eye    = VECTOR3(-4.30024, 1.52698, 3.03218);
    _lookAt = VECTOR3(-3.5036, 1.45884, 2.43159);
    _up     = VECTOR3(0.0632219, 0.997567, -0.029319);

    _pauseFrame = 200;

    return true;
  }
};

}

#endif
