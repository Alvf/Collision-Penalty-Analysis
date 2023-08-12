#ifndef SHELL_V_REGRESSION_H
#define SHELL_V_REGRESSION_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/TRIANGLE_MESH_FASTER.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Hyperelastic/Shell/STVK.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/BW_STRETCH.h"
#include "Hyperelastic/Shell/BW_SHEAR.h"
#include "Hyperelastic/Shell/BARAFF_WITKIN.h"
#include "Hyperelastic/Shell/BENDING_SPRING.h"
#include "Hyperelastic/Shell/QUADRATIC_F_BENDING.h"
#include "Hyperelastic/Shell/DIHEDRAL.h"
#include "Collision/C_PLANES_EE.h"
#include <filesystem>

namespace HOBAK {

class SHELL_V_REGRESSION : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Shell V bending scene, regression testing" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;

    _triangleMeshFilename = string("../data/objs/v_cloth.obj");

    TRIANGLE_MESH::readObjFile(_triangleMeshFilename.c_str(), vertices, triangles);
    _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);
    _gravity = VECTOR3(0, -1.0, 0.0);

    const REAL stretchingStiffness = 8.0;
    const REAL bendingStiffness = 1e5;
    _strechingEnergy = new SHELL::ARAP(stretchingStiffness, 0.0);
    _bendingEnergy = new SHELL::QUADRATIC_F_BENDING(bendingStiffness);

    // this will determine the MOV and JSON filenames
    char buffer[100];
    sprintf(buffer, "_S%.1f_B%.1f", stretchingStiffness, bendingStiffness);
    _sceneName = string("shell_bending_V_regression");

    _spacing = 1.0/(21.0 - 1.0);
    setCollisions(2000.0, _spacing / 8.0);
    _shellSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy, *_vfGeneral, *_eeGeneral);

    const REAL cubeScale = 8;
    VECTOR3 center0(-0.5 * cubeScale, -0.4 * cubeScale, 0.4 * cubeScale);
    _kinematicShapes.push_back(new CUBE(center0, cubeScale));
    _shellSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _shellSolver->addKinematicCollisionObject(_kinematicShapes[0]);
    
    _eye    = VECTOR3(3.151, 0.689313, 1.26046);
    _lookAt = VECTOR3(2.30528, 0.345273, 0.852525);
    _up     = VECTOR3(-0.245811, 0.929652, -0.274441);

    _worldCenter = _triangleMesh->getRestTranslation();
    _pauseFrame = 50;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override 
  {
    _shellSolver->externalForces().setZero();
    _shellSolver->addGravity(_gravity);
    _shellSolver->solve(verbose);

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    drawTriangleMesh(*_triangleMesh, true);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    drawKinematicConstraints(_triangleMesh, _shellSolver);

    drawVertex(*_triangleMesh, _arrowCounter);
  };
#endif

protected:
  REAL _spacing;
};

}

#endif
