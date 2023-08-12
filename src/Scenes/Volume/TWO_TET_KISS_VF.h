#ifndef TWO_TET_KISS_VF_H
#define TWO_TET_KISS_VF_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

class TWO_TET_KISS_VF : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " One tet is just barely kissing another in an vertex-edge collision." << endl;
    cout << " The integrator should just barely push them apart according to the" << endl;
    cout << " collision epsilon" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "two_tet_kiss_vf";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/two_tets_vertex_face_kiss.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, 0, 0);

    REAL E = 10.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH(vertices, tets);
    _hyperelastic = new VOLUME::ARAP(mu, lambda);
    
    // setup collision objects
    const REAL collisionEps = 0.03;
    setCollisions(mu * 10.0, collisionEps);

    // build the time integrator
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);
    _volumeSolver->setDt(1.0 / 30.0);

    VECTOR3 center0(0.0,-4.90, 0.0);
    VECTOR3 center1(0.0,6.9, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _kinematicShapes.push_back(new CUBE(center1, 10.0));
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[1]);
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    // set damping constants
    _volumeSolver->setRayeligh(0,0);
    const REAL dampingFactor = 0.01;
    _volumeSolver->collisionDampingBeta() = dampingFactor;

    _volumeSolver->vertexFaceSelfCollisionsOn() = true;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = false;

    _eye    = VECTOR3(-2.13557, 0.947023, 2.39579);
    _lookAt = VECTOR3(-1.35707, 0.964614, 1.76839);
    _up     = VECTOR3(-0.00368317, 0.999718, 0.0234604);

    _pauseFrame = 100;

    return true;
  }
  
};

}

#endif
