#ifndef STRAND_COLLISION_SCENE_H
#define STRAND_COLLISION_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"
#include "Collision/C_PLANES_EE.h"

namespace HOBAK {

class STRAND_COLLISION_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand collision scene" << endl;
    cout << "=====================================================================" << endl;
  }

  static void printMatlab(const VECTOR& p, const string name)
  {
    // print out some Matlab-ready output
    cout << " " << name.c_str() << " = [";
    for (unsigned int x = 0; x < p.size(); x++)
      cout << p[x] << " ";
    cout << "];" << endl;
  }

  static void printMatlab(const vector<VECTOR3>& p, const string name)
  {
    // print out some Matlab-ready output
    cout << " " << name.c_str() << " = [";
    for (unsigned int y = 0; y < 3; y++)
    {
      for (unsigned int x = 0; x < p.size(); x++)
        cout << p[x][y] << " ";
      if (y != 2)
        cout << ";" << endl << "\t";
      else
        cout << "];" << endl;
    }
  }

  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "strand_collision";

    //const int totalTopPoints = 5;
    //const int totalBottomPoints = 8;
    //const int totalTopPoints = 100;
    //const int totalTopPoints = 50;
    //const int totalTopPoints = 30;
    const int totalTopPoints = 20;
    //const int totalTopPoints = 25;
    const int totalBottomPoints = 25;
    const REAL deltaTop = 3.0 / (REAL)(totalTopPoints - 1);
    const REAL deltaBottom = 3.0 / (REAL)(totalBottomPoints - 1);
    assert(totalTopPoints >= 4);

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;

    for (int x = 0; x < totalTopPoints; x++)
      restVertices.push_back(VECTOR3(0,0, deltaTop * x));

    // vector of all the strands
    vector<vector<int> > strandIndices;

    // store the indices of individual strands
    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strandIndices.push_back(strand0);

    const int nextStart = restVertices.size();
    for (int x = 0; x < totalBottomPoints; x++)
      restVertices.push_back(VECTOR3(-1.5 + x * deltaBottom,1,1.5));
    
    //restVertices.push_back(VECTOR3(-0.5,1,1.5));
    //restVertices.push_back(VECTOR3(0.5,1,1.5));
    //restVertices.push_back(VECTOR3(1.5,1,1.5));

    // store the indices of individual strands
    vector<int> strand1;
    for (unsigned int x = nextStart ; x < restVertices.size(); x++)
      strand1.push_back(x);
    strandIndices.push_back(strand1);

    //std::vector<VECTOR3> vertices = restVertices;

    _gravity = VECTOR3(0,981,0);
    //_gravity = VECTOR3(0,0,0);
    _gravityPreloadingStep = 0;
    _totalGravityPreloadingSteps = 0;

    const REAL E = 1e7; // Taz settings
    //const REAL E = 3.727e10;
    //const REAL E = 3.727e8;
    //const REAL E = 3.727e6;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;

    #if 0 // Set 1 to see debug loader
    if (!readSOBJFile("../data/debug_points.sobj", restVertices, strandIndices)) { std::exit(1); }
    #endif
    _strandMesh = new STRAND_MESH(restVertices, strandIndices,
                                      E, nu, density, radiusA, radiusB);

    // create the integrator
#if 0
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
#else
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(100.0);
#endif
    vector<REAL> springArgs;
    springArgs.push_back(2025); //stiffness
    springArgs.push_back(0.09); //eps
    SIGNED_LEN_PLANES* lengthFuncEE = new SIGNED_LEN_PLANES();
    ENERGY_1D* normalSpringEE = new ENERGY_1D(springArgs);
    C_PLANES_EE* cFuncEE = new C_PLANES_EE();

    ENERGY_12D* edgeEdgeCollisions = new ENERGY_12D(normalSpringEE, cFuncEE, lengthFuncEE);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *edgeEdgeCollisions);
    _strandSolver->setDt(1.0 / 30.0);

    // constrain the first vertex
    VECTOR3 center0(0,0,0);
    //const REAL cubeScale = 0.5;
    const REAL cubeScale = 2.1;
    _kinematicShapes.push_back(new CUBE(center0, cubeScale * deltaTop));
    VECTOR3 center1(0,0,3);
    _kinematicShapes.push_back(new CUBE(center1, cubeScale * deltaTop));
    VECTOR3 center2(-1.5,1,1.5);
    _kinematicShapes.push_back(new CUBE(center2, cubeScale * deltaBottom));
    VECTOR3 center3(1.5,1,1.5);
    _kinematicShapes.push_back(new CUBE(center3, cubeScale * deltaBottom));

    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    _eye    = VECTOR3(2.01296, -0.103389, -2.35256);
    _lookAt = VECTOR3(1.57852, 0.074342, -1.46956);
    _up     = VECTOR3(-0.0702835, -0.984038, 0.16349);

    _eye    = VECTOR3(2.3278, 2.05146, -2.63138);
    _lookAt = VECTOR3(1.89337, 2.22919, -1.74839);
    _up     = VECTOR3(-0.0702828, -0.984039, 0.163489);

    _eye    = VECTOR3(2.63797, 1.03746, -3.00624);
    _lookAt = VECTOR3(2.17208, 1.21477, -2.13935);
    _up     = VECTOR3(-0.0702833, -0.984038, 0.16349);

    _eye    = VECTOR3(2.30806, -0.305108, -2.46584);
    _lookAt = VECTOR3(1.84931, -0.102653, -1.60064);
    _up     = VECTOR3(-0.0820219, -0.97919, 0.185639);

    _pauseFrame = 460;
    //_pauseFrame = 315;
    //_pauseFrame = 91;

    return true;
  }

  void preload(const bool verbose)
  {
    cout << " Gravity preloading step " << _gravityPreloadingStep << " of " << _totalGravityPreloadingSteps << endl;
    const REAL gravityFraction = (float)_gravityPreloadingStep / (float) _totalGravityPreloadingSteps;
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(gravityFraction * _gravity);
    _strandSolver->solveDynamics(verbose);

    _gravityPreloadingStep++;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    const VECTOR3 translation(0, -0.01, 0.0); // 30 explodes
    //const VECTOR3 translation(0, -0.02, 0.0); // 20 explodes

    _kinematicShapes[0]->translation() -= translation;
    _kinematicShapes[1]->translation() -= translation;

    cout << " Translation: " << _kinematicShapes[0]->translation().transpose() << endl;

    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    //_strandSolver->solveDynamics(verbose);
    _strandSolver->solveNewton(verbose);

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    //drawAxes();
    drawStrandMesh(*_strandMesh);
    //drawStrandMeshOld(*_strandMesh);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    //drawCollisionsOld(*_strandMesh);
    drawCollisions(*_strandMesh);

    /*
    drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_triangleMesh, _strandSolver);
    drawStrandTwistFreeFrames(*_strandMesh);
    drawStrandTwistFrames(*_strandMesh);
    drawStrandTwistForces(*_strandMesh);
    */
  };
#endif

protected:
  //STRAND_MESH* _strandMesh;
  STRAND_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  int _gravityPreloadingStep;
  int _totalGravityPreloadingSteps;
};

}

#endif
