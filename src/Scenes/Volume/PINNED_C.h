#ifndef PINNED_C_H
#define PINNED_C_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

class PINNED_C : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " An C-shaped mesh is pinned to the ground, and when its top lip" << endl;
    cout << " collides with its bottom lip, a smooth self-collision response" << endl;
    cout << " should occur " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "pinned_c";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/c_10.tobj");
    //string testFile("../data/c_20.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

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
    _hyperelastic = new VOLUME::SNH(mu, lambda);

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
    
    // setup collision objects
    setCollisions(1000.0, 0.01);

    // build the time integrator
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);
    _volumeSolver->setDt(1.0 / 60.0);

    VECTOR3 center0(0.0,-4.77, 0.0);
    VECTOR3 center1(0.0,6.9, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _kinematicShapes.push_back(new CUBE(center1, 10.0));
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    // collision constants
    //const REAL collisionMu = 1000.0;  // default
    //_volumeSolver->collisionStiffness() = collisionMu;
    _volumeSolver->collisionDampingBeta() = 0.01;   //default

    _volumeSolver->vertexFaceSelfCollisionsOn() = true;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = true;

    // front view, see the two lips meeting
    _eye    = VECTOR3(-0.517937, 0.517941, 0.452489);
    _lookAt = VECTOR3(0.478759, 0.451224, 0.498817);
    _up     = VECTOR3(0.0662986, 0.997745, 0.0105329);

    // side view, see the jitter when the lips slide
    _eye    = VECTOR3(0.320933, 0.559737, -0.261081);
    _lookAt = VECTOR3(0.248036, 0.516995, 0.735343);
    _up     = VECTOR3(0.0265707, 0.998644, 0.044781);

    _worldCenter = VECTOR3(0.497, 0.785889, 0.452556);
    _pauseFrame = 400;

    return true;
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    SIMULATION_SCENE::drawScene();

    if (_drawFeature)
    {
      glDisable(GL_DEPTH_TEST);
      glPointSize(10.0);
      drawVertexFacePairs(*_tetMesh, _arrowCounter);
      //drawVertexFacePairs(*_tetMesh);
    }
  };
#endif
};

}

#endif
