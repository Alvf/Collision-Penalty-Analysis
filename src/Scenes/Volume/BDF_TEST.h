#ifndef BDF_TEST_H
#define BDF_TEST_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Timestepper/Volume/BDF_1.h"
#include "Timestepper/Volume/BDF_2.h"

namespace HOBAK {

class BDF_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Dropping a cylinder on the kinematic ground, making sure that it" << endl;
    cout << " reboundes correctly " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "bdf_test";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/cylinder_10.tobj");

    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, -1.0, 0);

    const REAL E = 3.0;  // default
    const REAL nu = 0.48; // lambda \approx 25  // default
    const REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    const REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(mu, lambda);

    setCollisions(1000.0, 0.01);

    // build the time integrator
    //_volumeSolver = new VOLUME::BDF_1(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);
    _volumeSolver = new VOLUME::BDF_2(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);

    _volumeSolver->vertexFaceSelfCollisionsOn() = false;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = false;

    const VECTOR3 center0(0.0,-5.0, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    _eye    = VECTOR3(-2.5996, 0.52361, 0.286395);
    _lookAt = VECTOR3(-1.60313, 0.44686, 0.32046);
    _up     = VECTOR3(0.0765102, 0.997036, 0.00830762);

    // give it a big ballistic velocity
    VECTOR velocity(vertices.size() * 3);
    velocity.setZero();
    for (int x = 0; x < velocity.size() / 3; x++)
      velocity[3 * x + 1] = -3.0;

    _volumeSolver->velocity() = velocity;

    // if it's BDF-1, need to set the old velocity too
    VOLUME::BDF_1* setOld = dynamic_cast<VOLUME::BDF_1*>(_volumeSolver);
    if (setOld != NULL)
      setOld->velocityOld() = velocity;

    // if it's BDF-2, need to set the old velocity too
    VOLUME::BDF_2* setOlder = dynamic_cast<VOLUME::BDF_2*>(_volumeSolver);
    if (setOlder != NULL)
    {
      setOlder->velocityOld() = velocity;
      setOlder->velocityOlder() = velocity;
    }

    _pauseFrame = 250;

    return true;
  }
 
// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    glDisable(GL_DEPTH_TEST);
    glPointSize(10.0);
    drawVertices(*_tetMesh, _volumeSolver->constrainedNodes());
    drawPlaneConstraints(_tetMesh, _volumeSolver);
  };
#endif

};

}

#endif
